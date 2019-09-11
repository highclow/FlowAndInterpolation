// This is a wrapper for Ce Liu's Coarse2Fine optical flow implementation.
// It converts the contiguous image array to the format needed by the optical
// flow code. Handling conversion in the wrapper makes the cythonization
// simpler.
// Author: Deepak Pathak (c) 2016

#include "SplatMotionsWrapper.h"
#include "Image.h"
#include "SplatMotions.h"
using namespace std;

void SplatMotionsWrapper(double *vx, double *vy,
                         const double *vxForward, const double *vyForward,
                         const double *vxBackward, const double *vyBackward,
                         const double *Im1, const double *Im2,
                         int colType, int h, int w, int c, double t) {
  DImage ImFormatted1, ImFormatted2;
  DImage vxForwardFormatted, vyForwardFormatted;
  DImage vxBackwardFormatted, vyBackwardFormatted;
  DImage vxFormatted, vyFormatted;

  // format input in the format needed by backend
  vxForwardFormatted.allocate(w, h, 1);
  vyForwardFormatted.allocate(w, h, 1);
  vxBackwardFormatted.allocate(w, h, 1);
  vyBackwardFormatted.allocate(w, h, 1);
  ImFormatted1.allocate(w, h, c);
  ImFormatted2.allocate(w, h, c);

  memcpy(vxForwardFormatted.pData, vxForward, h * w * sizeof(double));
  memcpy(vyForwardFormatted.pData, vyForward, h * w * sizeof(double));
  memcpy(vxBackwardFormatted.pData, vxBackward, h * w * sizeof(double));
  memcpy(vyBackwardFormatted.pData, vyBackward, h * w * sizeof(double));
  memcpy(ImFormatted1.pData, Im1, h * w * c * sizeof(double));
  memcpy(ImFormatted2.pData, Im2, h * w * c * sizeof(double));

  vxForwardFormatted.setColorType(1);
  vyForwardFormatted.setColorType(1);
  vxBackwardFormatted.setColorType(1);
  vyBackwardFormatted.setColorType(1);
  ImFormatted1.setColorType(colType);
  ImFormatted2.setColorType(colType);

  vxFormatted.allocate(w, h, 1);
  vyFormatted.allocate(w, h, 1);

  // call bidirection splat motions backend
  SplatMotions::splatMotionsBidirect(vxFormatted, vyFormatted,
                                     vxForwardFormatted, vyForwardFormatted,
                                     vxBackwardFormatted, vyBackwardFormatted,
                                     ImFormatted1, ImFormatted2, t);
//
//  // copy formatted output to a contiguous memory to be returned
  memcpy(vx, vxFormatted.pData, h * w * sizeof(double));
  memcpy(vy, vyFormatted.pData, h * w * sizeof(double));

  // clear c memory
  ImFormatted1.clear();
  ImFormatted2.clear();
  vxForwardFormatted.clear();
  vyForwardFormatted.clear();
  vxBackwardFormatted.clear();
  vyBackwardFormatted.clear();
  vxFormatted.clear();
  vyFormatted.clear();

  return;
}
