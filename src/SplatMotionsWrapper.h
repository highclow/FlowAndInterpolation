// This is a C++ reimplementation wrapper for booster library.
// Author: Lu He (c) 2019

// override-include-guard
extern void SplatMotionsWrapper(double *vx, double *vy,
                                const double *vxForward, const double *vyForward,
                                const double *vxBackward, const double *vyBackward,
                                const double * Im1, const double * Im2,
                                int colType, int h, int w, int c, int t);
