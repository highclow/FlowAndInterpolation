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

def image_to_cartesian(f, x, y):
    """ Converts coordinates in image space (zero at top left,
    reversed y) to cartesian space (with bounds at -0.5, 0.5) """
    h = f.shape[0]
    w = f.shape[1]
    ix = x / w - 0.5
    iy = 0.5 - y / h
    return (ix, iy)
    
def cartesian_to_image(f, x, y):
    """ Converts the cartesian coordinate(x, y) to image 
    coordinate(x, y). Returns actual coordinates -- NOT
    indices. """
    h = f.shape[0]
    w = f.shape[1]
    cx = (x + 0.5) * w
    cy = (0.5 - y) * h
    return (cx, cy)
    
def vector_cartesian_to_image(f, dx, dy):
    """ scales a vector (dx,dy) to image space dimensions. """
    h = f.shape[0]
    w = f.shape[1]
    return (w*dx, h*dy)
    
def vector_mage_to_cartesian(f, dx, dy):
    """ Scales a vector (dx, dy) to cartesian space dimensions. """
    h = f.shape[0]
    w = f.shape[1]
    return (dx/w, dy/h)
    
def cartesian_to_index(f, x, y):
    (ix, iy) = cartesian_to_image(f, x, y)
    return (int(round(ix)), int(round(iy)))

    
def splat(f, x, y):
    """ 
    Return the indices of pixels (as (x, y), not (row, col))
    of all pixels in a square 4-pixel splat centered at (x,y)
    """    
    h = f.shape[0]
    w = f.shape[1]
    center = cartesian_to_index(f, x, y)
    pixels = [[center[0], center[1]], [center[0], center[1]], [center[0], center[1]], [center[0], center[1]]]
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


def splat_v2(f, x, y):
    """ 
    Return the indices of pixels (as (x, y), not (row, col))
    of all pixels in a square 4-pixel splat centered at (x,y)
    """
    h = f.shape[0]
    w = f.shape[1]
    center = [int(round(x)), int(round(y))]
    pixels = [center for _ in range(4)]
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

    u should already be scaled.

    frame0 and frame1  should be in BGR.
    """
    h = frame0.shape[0]
    w = frame0.shape[1]
    ux = u[0]
    uy = u[1]
    xp0 = col - t*ux
    yp0 = row - t*uy
    xp1 = col + (1-t)*ux
    yp1 = row + (1-t)*uy

    xi0 = int(round(xp0))
    yi0 = int(round(yp0))
    xi1 = int(round(xp1))
    yi1 = int(round(yp1))
    if (check_indices(frame0, xi0, yi0) and
        check_indices(frame1, xi1, yi1)):
        i0=frame0[yi0][xi0]
        i1=frame1[yi1][xi1]
        return ut.eucdist(i0, i1)
    else:
        return 2 # default value if going OOB
        

def splat_forward(forward, frame0, frame1, splatty, pts, t=0.5):
    """ Splat the foward frame0 motion at [row][col] onto splatty. """
    h = frame0.shape[0]
    w = frame0.shape[1]
    # Scale to cartesian space
    motion_cart = forward.copy()
    motion_cart[:,:,0] = motion_cart[:,:,0] * t
    motion_cart[:,:,1] = motion_cart[:,:,1] * t
    for row in tqdm(range(h), position=True, desc="splatting forward"):
      for col in range(w):
        # xp and yp are the coordinates in the interpolated image
        xp = col + motion_cart[row][col][0]
        yp = row + motion_cart[row][col][1]
        splats = splat_v2(splatty, xp, yp)
        for s in splats:
          if check_indices(splatty, s[0], s[1]):
            motion = forward[row][col]
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
    motion_cart = backward.copy()
    motion_cart[:,:,0] = motion_cart[:,:,0] * (1-t)
    motion_cart[:,:,1] = motion_cart[:,:,1] * (1-t)
    for row in tqdm(range(h), position=True, desc="splatting backward"):
      for col in range(w):
        # Scale to cartesian space
        # xp and yp are the coordinates in the interpolated image.
        xp = col + motion_cart[row][col][0]
        yp = row + motion_cart[row][col][1]
        splats = splat_v2(splatty, xp, yp)
        for s in splats:
          if check_indices(splatty, s[0], s[1]):
            motion = -1*backward[row][col]
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


def check_bounds_cartesian(f, x, y):
    """ Returns true if (x,y) is in the cartesian bounds 
    of f, else returns false. """
    return (-0.5 <= x < 0.5) and (-0.5 <= y < 0.5)

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
    m0h = np.array([m0[0], m0[1], 1])
    m1h = np.array([m1[0], m1[1], 1])
    top = m0h.dot(m1h)
    
    return 1 - (top / (np.linalg.norm(m0h) * np.linalg.norm(m1h)))

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
