#ifndef __FLOW_INTERPOLATION_H__
#define __FLOW_INTERPOLATION_H__

#include "Image.h"

typedef double _FlowPrecision;

class FlowInterpolation
{
public:
    static void splatMotionsBidirect(DImage& vx, DImage& vy, const DImage& vxForward, const DImage& vyForward, const DImage& vxBackward, const DImage& vyBackward, const DImage& Im1, const DImage& Im2, double t);
    static void splatForward(DImage& vx, DImage& vy, DImage &pts, const DImage& vxForward, const DImage& vyForward, const DImage& Im1, const DImage& Im2, double t);
    static void splatBackward(DImage& vx, DImage& vy, DImage &pts, const DImage& vxBackward, const DImage& vyBackward, const DImage& Im1, const DImage& Im2, double t);
    static void findHoles(std::vector<int>& holes, DImage& vx, DImage& vy);
    static void fillHoles(std::vector<int>& holes, DImage& vx, DImage& vy);
    static void killMaxLimits(std::vector<int>& holes, DImage& vx, DImage& vy);
    static void colorTransfer(DImage& dest, const DImage& interp, const DImage& Im1, const DImage& Im2, const DImage& forward, const DImage& backward, double t);


    static inline void splat(std::vector<int> &pixels, double x, double y){
        int cx = static_cast<int>(x + 0.5);
        int cy = static_cast<int>(y + 0.5);
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

    static inline bool insideBoundary(int h, int w, int row, int col, int neighbor) {
       if ((row == 0 && neighbor < -1) || (row == h-1 && neighbor > 1) || 
          (col == 0 && neighbor % w == -1) || (col == w-1 && neighbor % w == 1)) {
         return false;
       }
       return true;
    }

    static inline bool checkIndices(int h, int w, int x, int y) {
      return 0 <= y and y < h and 0 <= x and x < w;
    }

    static inline double squareDistance(const double *&p1, const double *&p2, int len) {
      double sum = 0;
      double diff;
      for (int i=0; i!=len; ++i){
        diff = *p1 - *p2;
        sum += diff * diff;
        ++p1;
        ++p2;
      }
      return sum;
    }

    static inline double flowDistance(double fx, double fy, double ux, double uy) {
      double distance = fx * ux + fy * uy + 1;
      double norm = sqrt((fx*fx + fy*fy + 1) *  (ux*ux + uy*uy + 1));
      return 1 - distance / norm;
    }

    static inline double intensity(const DImage& Im1, const DImage& Im2, double ux, double uy, int row, int col, double t) {
      int h = Im1.height();
      int w = Im1.width();
      int c = Im1.nchannels();
      double xp0 = col - t * ux + 0.5;
      double yp0 = row - t * uy + 0.5;
      double xp1 = ux + xp0;
      double yp1 = uy + yp0;
      int xi0 = static_cast<int>(xp0);
      int yi0 = static_cast<int>(yp0);
      int xi1 = static_cast<int>(xp1);
      int yi1 = static_cast<int>(yp1);
      if (checkIndices(h, w, xi0, yi0) && checkIndices(h, w, xi1, yi1)) {
        const double *i0 = Im1.data();
        const double *i1 = Im2.data();
        i0 += (yi0 * w + xi0) * c;
        i1 += (yi1 * w + xi1) * c;
        return squareDistance(i0, i1, c);
      }
      return 4;
    }
};



#endif // __FLOW_INTERPOLATION_H__
