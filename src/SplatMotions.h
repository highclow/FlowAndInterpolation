#ifndef __SPLATMOTIONS_H__
#define __SPLATMOTIONS_H__

#include "Image.h"

typedef double _FlowPrecision;

class SplatMotions
{
public:
    static void splatMotionsBidirection(DImage& vx, DImage& vy, const DImage& vxForward, const DImage& vyForward, const DImage& vxBackward, const DImage& vyBackward, const DImage& Im1, const DImage& Im2, float t=0.5);


}



#endif // __SPLATMOTIONS_H__
