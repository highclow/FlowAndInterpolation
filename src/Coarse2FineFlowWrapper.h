// This is a wrapper for Ce Liu's Coarse2Fine optical flow implementation.
// It converts the contiguous image array to the format needed by the optical
// flow code. Handling conversion in the wrapper makes the cythonization
// simpler.
// Author: Deepak Pathak (c) 2016

// override-include-guard
extern void Coarse2FineFlowWrapper(double * vx, double * vy, double * warpI2,
                              const double * Im1, const double * Im2,
                              double alpha, double ratio, int minWidth,
                              int nOuterFPIterations, int nInnerFPIterations,
                              int nSORIterations, int colType,
                              int h, int w, int c);


extern void SplatMotionsWrapper(double *vx, double *vy,
                                const double *vxForward, const double *vyForward,
                                const double *vxBackward, const double *vyBackward,
                                const double * Im1, const double * Im2,
                                int colType, int h, int w, int c, double t);


extern void ColorTransferWrapper(double *dest, const double *flow,
                             const double *Im1, const double *Im2,
                             const double *forward, const double *backward,
                             int colType, int h, int w, int c, double t);

