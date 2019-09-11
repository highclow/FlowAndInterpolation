# Author: Deepak Pathak (c) 2016

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# from __future__ import unicode_literals
import numpy as np
from PIL import Image
import cv2
import time
import argparse
import os
import pyflow
import booster.pixels as pix
import booster.color_transfer as ct


# Flow Options:
alpha = 0.012
ratio = 0.5
minWidth = 20
nOuterFPIterations = 7
nInnerFPIterations = 1
nSORIterations = 30
colType = 0  # 0 or default:RGB, 1:GRAY (but pass gray image with shape (h,w,1))

def interpolation(im1, im2, ts):
  s = time.time()
  uForward, vForward, im2WForward = pyflow.coarse2fine_flow(im1, im2,
          alpha, ratio, minWidth, nOuterFPIterations, nInnerFPIterations,
          nSORIterations, colType)
  uBackward, vBackward, im2WBackward = pyflow.coarse2fine_flow(im2, im1,
          alpha, ratio, minWidth, nOuterFPIterations, nInnerFPIterations,
          nSORIterations, colType)
  e = time.time()
  print('Time Taken: %.2f seconds for image of size (%d, %d, %d)' % (
          e - s, im1.shape[0], im1.shape[1], im1.shape[2]))
  s = time.time()
  interpolated = []
  for t in ts:
    flow = pyflow.splat_motions(uForward, vForward, uBackward, vBackward,
                                im1, im2, t)
    flow = np.concatenate((flow[0][..., None], flow[1][..., None]), axis=2)
    forward = np.concatenate((uForward[..., None], vForward[..., None]), axis=2)
    backward = np.concatenate((uBackward[..., None], vBackward[..., None]), axis=2)
    pix.fill_holes(flow)
    pix.fill_holes(flow)
    pix.kill_infs(flow)
    flow = cv2.GaussianBlur(flow, (11, 11), 10)
    interp = pyflow.color_transfer(im1,
                                   im2,
                                   forward,
                                   backward,
                                   flow,
                                   t)
    interpolated.append(interp)
  e = time.time()
  print('Time Taken: %.2f seconds for image of size (%d, %d, %d)' % (
          e - s, im1.shape[0], im1.shape[1], im1.shape[2]))
  return interpolated

if __name__ == '__main__':
  im1 = cv2.imread('examples/w1SyWlV9444_00139_s000-02081.jpg')
  im2 = cv2.imread('examples/w1SyWlV9444_00139_s000-00411.jpg')
  im1 = im1.astype(float) / 255.
  im2 = im2.astype(float) / 255.
# Interpolate rate
  ts = np.arange(0.0, 1.0, 0.1)
  interpolated = interpolation(im1, im2, ts)
  for k, interp in zip(ts,interpolated):
    cv2.imwrite('examples/interpolated_%02d.jpg'%(int(10*k)), (interp * 255).astype(int))
