#include <limits>
#include "Matrix.h"
#include "FlowInterpolation.h"


void FlowInterpolation::fillHoles(DImage& vx, DImage& vy) {
    const double *pvx = vx.data();
    const double *pvy = vy.data();
    const int nPixels = vx.npixels();
    std::vector<int> indices(nPixels); 
    std::vector<double> sum(nPixels);
    std::vector<int> count(nPixels);
    int cnt = 0;
    for (int i=0; i != nPixels; ++i) {
      if (pvx[i] == std::numeric_limits<double>::max() || 
          pvy[i] == std::numeric_limits<double>::max()) {
        indices[cnt] = i;
//        for (int m1=-1; m1!=2; ++y)
      }
    }
}

void FlowInterpolation::splatForward(DImage& vx, DImage& vy, DImage &pts, const DImage& vxForward, const DImage& vyForward, const DImage& Im1, const DImage& Im2, double  t){
    
    const int nRows = vxForward.height();
    const int nCols = vxForward.width();
    const int nChannels = Im2.nchannels();

    const double *pvxForward = vxForward.data();
    const double *pvyForward = vyForward.data();
    const double *pIm1 = Im1.data();
    const double *pIm2 = Im2.data();

    double ux, uy, xp, yp;
    std::vector<int> pixels(8);
    for (int i=0; i != nRows; ++i) {
      for (int j=0; j != nCols; ++j) {
        ux = *pvxForward;
        uy = *pvyForward; 
        xp = j + t * ux;
        yp = i + t * uy;
        splat(pixels, xp, yp);

        for (int k=0; k != 4; ++k) {
          int col = k << 1;
          double px = pixels[col];
          double py = pixels[col+1];
          if (checkIndices(nRows, nCols, px, py)) {
            double new_ptc = intensity(Im1, Im2, ux, uy, py, px, t);
            int offset = py * nCols + px;
            if (new_ptc < pts[offset]){
              pts[offset] = new_ptc;
              vx[offset] = ux;
              vy[offset] = uy;
            }
          }
        }
        ++pvxForward;
        ++pvyForward;
        pIm1 += nChannels;
        pIm2 += nChannels;
      }
    }
}


void FlowInterpolation::splatBackward(DImage& vx, DImage& vy, DImage &pts, const DImage& vxBackward, const DImage& vyBackward, const DImage& Im1, const DImage& Im2, double  t){

    const int nRows = vxBackward.height();
    const int nCols = vxBackward.width();
    const int nChannels = Im2.nchannels();

    const double *pvxBackward = vxBackward.data();
    const double *pvyBackward = vyBackward.data();
    const double *pIm1 = Im1.data();
    const double *pIm2 = Im2.data();

    double ux, uy, xp, yp;
    std::vector<int> pixels(8);
    for (int i=0; i != nRows; ++i) {
      for (int j=0; j != nCols; ++j) {
        ux = *pvxBackward;
        uy = *pvyBackward;
        xp = j + (1-t) * ux;
        yp = i + (1-t) * uy;
        splat(pixels, xp, yp);

        for (int k=0; k != 4; ++k) {
          int col = k << 1;
          double px = pixels[col];
          double py = pixels[col+1];
          if (checkIndices(nRows, nCols, px, py)) {
            double new_ptc = intensity(Im1, Im2, -ux, -uy, py, px, t);
            int offset = py * nCols + px;
            if (new_ptc < pts[offset]){
              pts[offset] = new_ptc;
              vx[offset] = -ux;
              vy[offset] = -uy;
            }
          }
        }
        ++pvxBackward;
        ++pvyBackward;
        pIm1 += nChannels;
        pIm2 += nChannels;
      }
    }
}
    


void FlowInterpolation::splatMotionsBidirect(DImage& vx, DImage& vy, const DImage& vxForward, const DImage& vyForward, const DImage& vxBackward, const DImage& vyBackward, const DImage& Im1, const DImage& Im2, double t){

    const int nRows = vxForward.height();
    const int nCols = vxForward.width();
    DImage pts;
    pts.setValue(5.0, nRows, nCols, 1);
    splatForward(vx, vy, pts, vxForward, vyForward, Im1, Im2, t);
    splatBackward(vx, vy, pts, vxBackward, vyBackward, Im1, Im2, t);
    

//    # Fill holes twice to get rid of big holes
//    fill_holes(splatty)
//    fill_holes(splatty)
//    kill_nans(splatty)
  
}


void FlowInterpolation::colorTransfer(DImage& dest, const DImage& interp, const DImage& Im1, const DImage& Im2, const DImage& forward, const DImage& backward, double t) {

    const double *pInterp = interp.data();
    double *pDest = dest.data();
    const int nRows = Im1.height();
    const int nCols = Im1.width();
    const int nChannels = Im1.nchannels();
    double ux, uy, xp1, yp1, xp2, yp2;
    int xi1, yi1, xi2, yi2;

    for (int row=0; row!=nRows; ++row) {
      for (int col=0; col!=nCols; ++col) {
        ux = *pInterp;
        uy = *(pInterp+1);
        xp1 = col - t * ux + 0.5;
        yp1 = row - t * uy + 0.5;
        xp2 = ux + xp1;
        yp2 = uy + yp1;
        xi1 = static_cast<int>(xp1);
        yi1 = static_cast<int>(yp1);
        xi2 = static_cast<int>(xp2);
        yi2 = static_cast<int>(yp2);
        if (checkIndices(nRows, nCols, xi1, yi1) &&
            checkIndices(nRows, nCols, xi2, yi2)) {
          int paddr1 = yi1 * nCols + xi1;
          int paddr2 = yi2 * nCols + xi2;
          int offset1 = paddr1 << 1;
          int offset2 = paddr2 << 1;
          double fx = forward[offset1];
          double fy = forward[offset1+1];
          double bx = -backward[offset2];
          double by = -backward[offset2+1];
          double d1 = flowDistance(fx, fy, ux, uy);
          double d2 = flowDistance(bx, by, ux, uy);
          double diff = abs(d1 - d2);
          offset1 = paddr1 * 3;
          offset2 = paddr2 * 3;
          if (diff < 0.1) {
            pDest[0] = (Im1[offset1] + Im2[offset2]) * 0.5;
            pDest[1] = (Im1[offset1+1] + Im2[offset2+1]) * 0.5;
            pDest[2] = (Im1[offset1+2] + Im2[offset2+2]) * 0.5;
          } else if (d1 > d2) {
            pDest[0] = Im2[offset2];
            pDest[1] = Im2[offset2+1];
            pDest[2] = Im2[offset2+2];
          } else {
            pDest[0] = Im1[offset1];
            pDest[1] = Im1[offset1+1];
            pDest[2] = Im1[offset1+2];
          }
        }
        pInterp += 2;
        pDest += 3;
      }
    }
    return;
}
