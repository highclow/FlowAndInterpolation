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
minWidth = 40
nOuterFPIterations = 7
nInnerFPIterations = 1
nSORIterations = 30
colType = 0  # 0 or default:RGB, 1:GRAY (but pass gray image with shape (h,w,1))
# Interpolate rate
t = 0.5 #

def interpolation(im1, im2):
  uForward, vForward, im2WForward = pyflow.coarse2fine_flow(im1, im2,
          alpha, ratio, minWidth, nOuterFPIterations, nInnerFPIterations,
          nSORIterations, colType)
  uBackward, vBackward, im2WBackward = pyflow.coarse2fine_flow(im2, im1,
          alpha, ratio, minWidth, nOuterFPIterations, nInnerFPIterations,
          nSORIterations, colType)
  forward = np.concatenate((uForward[..., None], vForward[..., None]), axis=2)
  backward = np.concatenate((uBackward[..., None], vBackward[..., None]), axis=2)
  flow = pix.splat_motions_bidi(forward, backward,
                                im1, im2, t)
  flow = cv2.GaussianBlur(flow, (11, 11), 10)
  interpolated = ct.color_transfer_occlusions(im1,
                                              im2,
                                              forward,
                                              backward,
                                              flow,
                                              t)
  return interpolated

if __name__ == '__main__':
  im1 = cv2.imread('examples/car1.jpg')
  im2 = cv2.imread('examples/car2.jpg')
  im1 = im1.astype(float) / 255.
  im2 = im2.astype(float) / 255.
  s = time.time()
  interpolated = interpolation(im1, im2)
  e = time.time()
  print('Time Taken: %.2f seconds for image of size (%d, %d, %d)' % (
          e - s, im1.shape[0], im1.shape[1], im1.shape[2]))
  cv2.imwrite('examples/interpolated.jpg', (interpolated * 255).astype(int))
