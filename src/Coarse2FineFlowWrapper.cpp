// This is a wrapper for Ce Liu's Coarse2Fine optical flow implementation.
// It converts the contiguous image array to the format needed by the optical
// flow code. Handling conversion in the wrapper makes the cythonization
// simpler.
// Author: Deepak Pathak (c) 2016

#include <limits>
#include "Coarse2FineFlowWrapper.h"
#include "Image.h"
#include "OpticalFlow.h"
#include "FlowInterpolation.h"
using namespace std;


void Coarse2FineFlowWrapper(double * vx, double * vy, double * warpI2,
                              const double * Im1, const double * Im2,
                              double alpha, double ratio, int minWidth,
                              int nOuterFPIterations, int nInnerFPIterations,
                              int nSORIterations, int colType,
                              int h, int w, int c) {
  DImage ImFormatted1, ImFormatted2;
  DImage vxFormatted, vyFormatted, warpI2Formatted;

  // format input in the format needed by backend
  ImFormatted1.allocate(w, h, c);
  ImFormatted2.allocate(w, h, c);
  memcpy(ImFormatted1.pData, Im1, h * w * c * sizeof(double));
  memcpy(ImFormatted2.pData, Im2, h * w * c * sizeof(double));
  ImFormatted1.setColorType(colType);
  ImFormatted2.setColorType(colType);

  // call optical flow backend
  OpticalFlow::Coarse2FineFlow(vxFormatted, vyFormatted, warpI2Formatted,
                                ImFormatted1, ImFormatted2,
                                alpha, ratio, minWidth,
                                nOuterFPIterations, nInnerFPIterations,
                                nSORIterations);

  // copy formatted output to a contiguous memory to be returned
  memcpy(vx, vxFormatted.pData, h * w * sizeof(double));
  memcpy(vy, vyFormatted.pData, h * w * sizeof(double));
  memcpy(warpI2, warpI2Formatted.pData, h * w * c * sizeof(double));

  // clear c memory
  ImFormatted1.clear();
  ImFormatted2.clear();
  vxFormatted.clear();
  vyFormatted.clear();
  warpI2Formatted.clear();

  return;
}


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
  vxFormatted.setValue(std::numeric_limits<double>::max(), w, h);
  vyFormatted.setValue(std::numeric_limits<double>::max(), w, h);

  ImFormatted1.setColorType(colType);
  ImFormatted2.setColorType(colType);


  // call bidirection splat motions backend
  FlowInterpolation::splatMotionsBidirect(vxFormatted, vyFormatted,
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


void ColorTransferWrapper(double *dest, const double *flow,
                          const double *Im1, const double *Im2,
                          const double *forward, const double *backward,
                          int colType, int h, int w, int c, double t) {
  
  DImage flowFormatted, destFormatted;
  DImage ImFormatted1, ImFormatted2;
  DImage forwardFormatted, backwardFormatted;

  // format input in the format needed by backend
  destFormatted.allocate(w, h, 3);
  flowFormatted.allocate(w, h, 2);
  ImFormatted1.allocate(w, h, c);
  ImFormatted2.allocate(w, h, c);
  forwardFormatted.allocate(w, h, 2);
  backwardFormatted.allocate(w, h, 2);

  memcpy(destFormatted.pData, dest, h * w * 3 * sizeof(double));
  memcpy(flowFormatted.pData, flow, h * w * 2 * sizeof(double));
  memcpy(ImFormatted1.pData, Im1, h * w * c * sizeof(double));
  memcpy(ImFormatted2.pData, Im2, h * w * c * sizeof(double));
  memcpy(forwardFormatted.pData, forward, h * w * 2 * sizeof(double));
  memcpy(backwardFormatted.pData, backward, h * w * 2 * sizeof(double));

  destFormatted.setColorType(colType);
  ImFormatted1.setColorType(colType);
  ImFormatted2.setColorType(colType);

  // call bidirection splat motions backend
  FlowInterpolation::colorTransfer(destFormatted, flowFormatted,
                                   ImFormatted1, ImFormatted2,
                                   forwardFormatted, backwardFormatted, t);
//
//  // copy formatted output to a contiguous memory to be returned
  memcpy(dest, destFormatted.pData, h * w * 2 * sizeof(double));

  // clear c memory
  destFormatted.clear();
  flowFormatted.clear();
  ImFormatted1.clear();
  ImFormatted2.clear();
  forwardFormatted.clear();
  backwardFormatted.clear();

  return;
}
