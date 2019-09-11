#include "SplatMotions.h"
#include "Matrix.h"


void SplatMotions::splatForward(DImage& vx, DImage& vy, DImage &pts, const DImage& vxForward, const DImage& vyForward, const DImage& Im1, const DImage& Im2, double  t){
    
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


void SplatMotions::splatBackward(DImage& vx, DImage& vy, DImage &pts, const DImage& vxBackward, const DImage& vyBackward, const DImage& Im1, const DImage& Im2, double  t){

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
    


void SplatMotions::splatMotionsBidirect(DImage& vx, DImage& vy, const DImage& vxForward, const DImage& vyForward, const DImage& vxBackward, const DImage& vyBackward, const DImage& Im1, const DImage& Im2, double t){

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
