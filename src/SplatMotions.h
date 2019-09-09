#ifndef __SPLATMOTIONS_H__
#define __SPLATMOTIONS_H__

#include "Image.h"

typedef double _FlowPrecision;

class SplatMotions
{
public:
    static void splatMotionsBidirect(DImage& vx, DImage& vy, const DImage& vxForward, const DImage& vyForward, const DImage& vxBackward, const DImage& vyBackward, const DImage& Im1, const DImage& Im2, float t=0.5);
    static void splatForward(DImage& vx, DImage& vy, UCImage &flag, DImage &pts, const DImage& vxForward, const DImage& vyForward, const DImage& Im1, const DImage& Im2, float t=0.5);


    static inline void splat(std::vector<int> &pixels, int h, int w, int x, int y){
        int cx = static_cast<int>(round((x + 0.5) * w));
        int cy = static_cast<int>(round((-y + 0.5) * h));
        pixels[0] = pixels[4] = cx;
        if (x > cx) {
            pixels[2] = pixels[6] = cx + 1;
        } else {
            pixels[2] = pixels[6] = cx - 1;
        }
        pixels[1] = pixels[3] = cy;
        if (y > cy) {
            pixels[5] = pixels[7] = cy + 1;
        } else {
            pixels[5] = pixels[7] = cy - 1;
        }
    }

    static inline bool checkIndices(int h, int w, int x, int y) {
      return 0 <= y and y < h and 0 <= x and x < w;
    }
};



#endif // __SPLATMOTIONS_H__
