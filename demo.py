# Author: Deepak Pathak (c) 2016

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# from __future__ import unicode_literals
import numpy as np
import cv2
import time
import argparse
import os
import pyflow

parser = argparse.ArgumentParser(
    description='Demo for python wrapper of Coarse2Fine Optical Flow')
parser.add_argument(
    '-input', dest='inv', default='/data/lu/open_dataset/FaceAntiSpoofing/SiW/SiW_release/Test/live/001/001-1-1-1-1.mov',
    #'-input', dest='inv', default='/data/lu/open_dataset/FaceAntiSpoofing/SiW/SiW_release/Test/spoof/001/001-1-2-1-1.mov',
    help='Input Videos')
args = parser.parse_args()

#im1 = np.array(Image.open('live/00311110450.jpg'))
#im2 = np.array(Image.open('live/00311110451.jpg'))
#im1 = np.array(Image.open('spoof/00312110169.jpg'))
#im2 = np.array(Image.open('spoof/00312110170.jpg'))
#im1 = np.array(Image.open('spoof/00323310393.jpg'))
#im2 = np.array(Image.open('spoof/00323310394.jpg'))
#im1 = np.array(Image.open('spoof/00313310393.jpg'))
#im2 = np.array(Image.open('spoof/00313310394.jpg'))
#im1 = np.array(Image.open('spoof/00313110168.jpg'))
#im2 = np.array(Image.open('spoof/00313110169.jpg'))

# Flow Options:
alpha = 0.012
ratio = 0.7
minWidth = 20
nOuterFPIterations = 7
nInnerFPIterations = 1
nSORIterations = 30
colType = 0  # 0 or default:RGB, 1:GRAY (but pass gray image with shape (h,w,1))

#cap = cv2.VideoCapture('/data/lu/open_dataset/FaceAntiSpoofing/SiW/SiW_release/Test/live/001/001-1-1-1-2.mov')
videoname = args.inv
annofile = videoname.replace("mov","face")
cap = cv2.VideoCapture(videoname)
save_prefix = args.inv.split('/')[-1].split('.')[0]
save_path = os.path.join('examples', '/'.join(args.inv.split('/')[-4:-1]))
if not os.path.exists(save_path):
  os.makedirs(save_path)
annotation = [list(map(int, d.split())) for d in open(annofile).readlines()]

def get_bbox(boxes):
    boxes = [box for box in boxes if box[2]!=0 and box[3]!=0]
    if len(boxes) == 0:
      return 0, 0, 0, 0
    x1, y1 = [min(box[i] for box in boxes) for i in range(2)]
    x2, y2 = [max(box[i] for box in boxes) for i in range(2,4)]
    height, width = y2 - y1, x2 - x1
    x1, y1 = int(x1 - width*0.5), int(y1 - height*0.5)
    x2, y2 = int(x2 + width*0.5), int(y2 + height*0.5)
    return max(0,x1//4), max(0,y1//4), x2//4, y2//4

ret, im1 = cap.read()
height, width, _ = im1.shape
im1 = cv2.resize(im1,(width//4, height//4))
im1 = im1.astype(float) / 255.
print(save_path, im1.shape)

cnt = 0
while True:
  ret, im2 = cap.read()
  if not ret:
    break
  im2 = cv2.resize(im2,(width//4, height//4))
  im2 = im2.astype(float) / 255.
  s = time.time()
  if cnt > 30 and cnt < 450 and (cnt+1) % 10  == 0:
    x1, y1, x2, y2 = get_bbox(annotation[cnt-2:cnt])
    if x2-x1 >= 40 and y2-y1 >= 40:
      u, v, im2W = pyflow.coarse2fine_flow(
        im1, im2,
        alpha, ratio, minWidth, nOuterFPIterations, nInnerFPIterations,
        nSORIterations, colType)
      e = time.time()
      print('Time Taken: %.2f seconds for image of size (%d, %d, %d)' % (
        e - s, im1.shape[0], im1.shape[1], im1.shape[2]))
      flow = np.concatenate((u[..., None], v[..., None]), axis=2)
      #np.save('examples/outFlow.npy', flow)
  
      hsv = np.zeros(im1.shape, dtype=np.uint8)
      hsv[:, :, 0] = 255
      hsv[:, :, 1] = 255
      mag, ang = cv2.cartToPolar(flow[..., 0], flow[..., 1])
      hsv[..., 0] = ang * 180 / np.pi / 2
      hsv[..., 2] = cv2.normalize(mag, None, 0, 255, cv2.NORM_MINMAX)
      rgb = cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR)
#      rgb = rgb[y1:y2,x1:x2].copy(order='C')
#      im2W = im2W[y1:y2,x1:x2].copy(order='C')
      cv2.imwrite('%s/%s-flow-%05d.jpg'%(save_path,save_prefix,cnt), rgb)
      cv2.imwrite('%s/%s-wrap-%05d.jpg'%(save_path,save_prefix,cnt), im2W * 255)

  im1 = im2
  cnt += 1
