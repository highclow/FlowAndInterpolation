"""
Functions for manipulating single pixels or patches of pixels
in an image.

"""
import time
import numpy as np
import cv2
import booster.utility as ut
from tqdm import tqdm

__author__     = "Henry Cooney"
__credits__    = ["Henry Cooney", "Feng Liu"]
__license__    = "MIT"
__version__    = "0.1"
__maintainer__ = "Henry Cooney"
__email__      = "hacoo36@gmail.com"
__status__     = "Prototype"
__repo__       = "https://github.com/hacoo/framebooster.git" 

    
def splat(f, x, y):
    """ 
    Return the indices of pixels (as (x, y), not (row, col))
    of all pixels in a square 4-pixel splat centered at (x,y)
    """
    h = f.shape[0]
    w = f.shape[1]
    center = [int(x+0.5), int(y+0.5)]
    pixels = [center.copy() for _ in range(4)]
    if x > center[0]:
        pixels[1][0] += 1
        pixels[3][0] += 1
    else:
        pixels[1][0] -= 1
        pixels[3][0] -= 1
    if y > center[1]:
        pixels[2][1] += 1
        pixels[3][1] += 1
    else:
        pixels[2][1] -= 1
        pixels[3][1] -= 1
    return pixels


def follow_intensity(frame0, frame1, u, 
                     row, col, t=0.5):
    """ Follow the flow u forward and backward from the location
    [row][col], and return the intensity difference between
    the forward and backward pixel. 

    frame0 and frame1  should be in BGR.
    """
    xp0 = col - t*u[0] + 0.5
    yp0 = row - t*u[1] + 0.5
    xp1 = u[0] + xp0
    yp1 = u[1] + yp0

    xi0 = int(xp0)
    yi0 = int(yp0)
    xi1 = int(xp1)
    yi1 = int(yp1)
    if (check_indices(frame0, xi0, yi0) and
        check_indices(frame1, xi1, yi1)):
        i0=frame0[yi0][xi0]
        i1=frame1[yi1][xi1]
        return ut.squaredist(i0, i1)
    else:
        return 4 # default value if going OOB
        

def splat_forward(forward, frame0, frame1, splatty, pts, t=0.5):
    """ Splat the foward frame0 motion at [row][col] onto splatty. """
    h = frame0.shape[0]
    w = frame0.shape[1]
    # Scale to cartesian space
    coords = forward.copy()
    coords[:,:,0] = coords[:,:,0] * t + np.arange(w)
    coords[:,:,1] = coords[:,:,1] * t + np.arange(h)[:,np.newaxis]
    coords = coords.reshape((h*w,2))
    motions = forward.reshape((h*w,2))
    for p, motion in tqdm(zip(coords, motions), position=True, desc="forward"):
        # xp and yp are the coordinates in the interpolated image
        splats = splat(splatty, p[0], p[1])
        for s in splats:
          if check_indices(splatty, s[0], s[1]):
            new_ptc = follow_intensity(frame0, frame1, motion, s[1], s[0], t)
            if new_ptc < pts[s[1]][s[0]]:
              splatty[s[1]][s[0]] = motion
              pts[s[1]][s[0]] = new_ptc


def splat_backward(backward, frame0, frame1, splatty, pts, t=0.5):
    """ 
    Splat the backward m motion at [row][col] onto splatty.
    frame0 should come first chronologically.
    
    """
    h = frame0.shape[0]
    w = frame0.shape[1]
    # Scale to cartesian space
    coords = backward.copy()
    coords[:,:,0] = coords[:,:,0] * (1-t) + np.arange(w)
    coords[:,:,1] = coords[:,:,1] * (1-t) + np.arange(h)[:,np.newaxis]
    coords = coords.reshape((h*w,2))
    motions = -backward.reshape((h*w,2))
    for p, motion in tqdm(zip(coords, motions), position=True, desc="backward"):
        # xp and yp are the coordinates in the interpolated image.
        splats = splat(splatty, p[0], p[1])
        for s in splats:
          if check_indices(splatty, s[0], s[1]):
            new_ptc = follow_intensity(frame0, frame1, motion, s[1], s[0], t)
            if new_ptc < pts[s[1]][s[0]]:
              splatty[s[1]][s[0]] = motion
              pts[s[1]][s[0]] = new_ptc


def splat_motions_bidi(forward, backward, frame0, frame1, t=0.5):
    """ Bidirectionally splat motions between frame 0 and frame1,
    with forward and backward flows. Return the interpolated flows. """
    h = frame0.shape[0]
    w = frame0.shape[1]
    splatty = np.zeros_like(forward)
    # Flows start out undefined
    splatty[:] = np.NAN
    ptcs = np.ones((h, w)) * float('inf')
    splat_forward(forward, frame0, frame1, splatty, ptcs, t)
    splat_backward(backward, frame0, frame1, splatty, ptcs, t)

    # Fill holes twice to get rid of big holes
    fill_holes(splatty)
    fill_holes(splatty)
    kill_nans(splatty)
    return splatty   
    

def kill_nans(f):
    """ Replaces all nans in the frame with [0,0] """
    h = f.shape[0]
    w = f.shape[1]
    for r in range(h):
        for c in range(w):
            if np.isnan(f[r][c][0]) or np.isnan(f[r][c][1]):
                f[r][c] = np.array([0.0, 0.0], dtype='float')

def kill_infs(f):
    """ Replaces all nans in the frame with [0,0] """
    h = f.shape[0]
    w = f.shape[1]
    for r in range(h):
        for c in range(w):
            if np.isinf(f[r][c][0]) or np.isinf(f[r][c][1]):
                f[r][c] = np.array([0.0, 0.0], dtype='float')
    
            
def fill_holes(splatty):
    """ Fill NaN holes in splatty by averaging neightbors. """
    indices = filter_not_nans(splatty)
    indices.sort(key=lambda x: num_nan_neighbors(splatty, x[0], x[1]))
    for i in tqdm(indices, position=True, desc="filling holes"): 
        average_fill(splatty, i)
    
        
def average_fill(f, indices):
    """ Fill in index i using the average of its neighbors. """
    n = 0
    u = 0.0
    v = 0.0
    for i in range(-1,2):
        for j in range(-1,2):
            r = indices[0] + i
            c = indices[1] + j
            if check_indices(f, c, r):
                pixel = f[r][c]
                if not np.isnan(pixel[0]) and not np.isnan(pixel[1]):
                    n += 1
                    u += pixel[0]
                    v += pixel[1]
    if n == 0:
        #averaged = np.array([0.0, 0.0], dtype='float')
        pass
    else:
        averaged = np.array([u/n, v/n], dtype='float')
        f[indices[0]][indices[1]] = averaged
    
def num_nan_neighbors(f, r, c):
    """ Return the number of neighbors of [r][c] that
    are NaN in the frame f. """
    nans = 0
    for i in range(-1,2):
        for j in range(-1,2):
            rp = r + i
            cp = c + j
            if check_indices(f, cp, rp):
                pixel = f[rp][cp]
                if np.isnan(pixel[0]) or np.isnan(pixel[1]):
                    nans += 1
    return nans
 

def filter_not_nans(f):
    """ Return a new list of indices, containing
    only pixels that are not nan in f. """
    h = f.shape[0]
    w = f.shape[1]
    no_nans = []
    for r in range(h):
        for c in range(w):
          pixel = f[r][c]
          if np.isnan(pixel[0]) or np.isnan(pixel[1]):
            no_nans.append((r,c))
    return no_nans


def check_indices(f, x, y):
    """ Check indices [y, x] """
    h = f.shape[0]
    w = f.shape[1]
    return (0 <= y < h) and (0 <= x < w)


def send_motion(src, row, col, t=0.5):
    """ find the pixel pointed to by the motion vector
    in src at [row, cow] """
    motion = src[row][col]
    return (int(round(row+t*motion[0])),
            int(round(col+t*motion[1])))
    
def send_motion_pixel(src, dest, row, col, t=0.5):
    """ Sends each motion vector in src to the corresponding
    closest pixel in dest, t timesteps ahead. """
    height = dest.shape[0]
    width  = dest.shape[1]
    dc = send_motion(src, row, col)
    if (dc[0] >= 0 and dc[0] < height) and (dc[1] >= 0 and dc[1] < width):
        dest[dc[0], dc[1]] = src[row][col]

def send_motions(src, t=0.5):
    """ Create a new with motion vectors from src, tranfered
    t timesteps into the future. """
    dest = np.zeros_like(src)
    rows = src.shape[0]
    cols = src.shape[1]
    for r in tqdm(range(rows)):
        for c in range(cols):
            send_motion_pixel(src, dest, r, c, t)
            return dest

def flowdist(m0, m1):
    """ Compute the distance between two optical flow vectors,
    m0 and m1. """

    # Convert to homogeneous coordinates, assuming timestep 
    # is just 1:
    assert(m0.shape == (2,))
    assert(m1.shape == (2,))
    return 1 - (np.dot(m0,m1)+1) / (np.sqrt((np.dot(m0,m0)+1) * (np.dot(m1,m1)+1)))

def patchdist(p0, p1):
    """ Compute the distance between two patches. Patches
    are 2d numpy arrays of flow vectors. """
    assert(p0.shape == p1.shape)
    D = p0.size
    sum = 0
    for i in range(p0.shape[0]):
        for j in range(p1.shape[1]):
            sum += flowdist(p0[i, j], p1[i, j])
    return sum / D

def makepatch(image, c, w):
    """ make a patch from the image, centered at c and with
    width w """
    rad = w // 2
    cx = c[0]
    cy = c[1]
    height = image.shape[0]
    width  = image.shape[1]
    patch  = np.zeros((w,w,image.shape[2]), dtype=image.dtype)
    
    for x in range(w):
            px = cx + x - rad
            for y in range(w):
                py = cy + y - rad
                if (0 <= py <= height) and (0 <= px <= width):
                    patch[y, x] = image[py, px]
    return patch
