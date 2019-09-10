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
# Interpolate rate
t = 0.7 #

def interpolation(im1, im2):
  s = time.time()
  uForward, vForward, im2WForward = pyflow.coarse2fine_flow(im1, im2,
          alpha, ratio, minWidth, nOuterFPIterations, nInnerFPIterations,
          nSORIterations, colType)
  uBackward, vBackward, im2WBackward = pyflow.coarse2fine_flow(im2, im1,
          alpha, ratio, minWidth, nOuterFPIterations, nInnerFPIterations,
          nSORIterations, colType)
  forward = np.concatenate((uForward[..., None], vForward[..., None]), axis=2)
  backward = np.concatenate((uBackward[..., None], vBackward[..., None]), axis=2)
  e = time.time()
  print('Time Taken: %.2f seconds for image of size (%d, %d, %d)' % (
          e - s, im1.shape[0], im1.shape[1], im1.shape[2]))
  s = time.time()
#  forward = np.load('./forward.npy')
#  backward = np.load('./backward.npy')
  flow = pix.splat_motions_bidi(forward, backward,
                                im1, im2, t)
#  uForward = forward[:,:,0].copy(order='C')
#  vForward = forward[:,:,1].copy(order='C')
#  uBackward = backward[:,:,0].copy(order='C')
#  vBackward = backward[:,:,1].copy(order='C')
#  flow2 = pyflow.splat_motions(uForward, vForward, uBackward, vBackward,
#                               im1, im2, t)
  flow = cv2.GaussianBlur(flow, (11, 11), 10)
  interpolated = ct.masked_transfer_color(im1,
                                          im2,
                                          forward,
                                          backward,
                                          flow,
                                          t)
  e = time.time()
  print('Time Taken: %.2f seconds for image of size (%d, %d, %d)' % (
          e - s, im1.shape[0], im1.shape[1], im1.shape[2]))
  return interpolated

if __name__ == '__main__':
  im1 = cv2.imread('examples/w1SyWlV9444_00139_s000-02081.jpg')
  im2 = cv2.imread('examples/w1SyWlV9444_00139_s000-00411.jpg')
  im1 = im1.astype(float) / 255.
  im2 = im2.astype(float) / 255.
  interpolated = interpolation(im1, im2)
  cv2.imwrite('examples/interpolated.jpg', (interpolated * 255).astype(int))
