#include <algorithm>
#include <limits>
#include "Matrix.h"
#include "FlowInterpolation.h"


std::ostream& operator<<(std::ostream& os, const std::vector<int> &input)
{
	for (auto const& i: input) {
		os << i << " ";
	}
	return os;
}



void FlowInterpolation::killMaxLimits(std::vector<int>& holes, DImage& vx, DImage& vy) {
    for (int i=0; i!=holes.size(); ++i) {
      vx[i] = 0.0;
      vy[i] = 0.0;
    }
}


void FlowInterpolation::findHoles(std::vector<int>& holes, DImage& vx, DImage& vy) {
    const double *pvx = vx.data();
    const double *pvy = vy.data();
    const int nPixels = vx.npixels();
    for (int k=0; k!=nPixels; ++k) {
      if (*pvx == std::numeric_limits<double>::max() ||
          *pvy == std::numeric_limits<double>::max()) {
        holes.push_back(k);
      }
      ++pvx;
      ++pvy;
    }
}


void FlowInterpolation::fillHoles(std::vector<int>& holes, DImage& vx, DImage& vy) {
    const int nRows = vx.height();
    const int nCols = vx.width();
    const int size = holes.size();
    int neighbors[8] = {-nCols-1, -nCols, -nCols+1, -1, 1, nCols-1, nCols, nCols+1};
    for (int i=0; i!=size; ++i) {
      int ind = holes[i];
      int row = ind / nCols;
      int col = ind - row * nCols;
      double sumx = 0.0;
      double sumy = 0.0;
      int cnt = 0;
      for (int k=0; k!=8; ++k) {
        if (insideBoundary(nRows, nCols, row, col, neighbors[k])) {
          int offset = ind + neighbors[k];
          if (vx[offset] != std::numeric_limits<double>::max() ||
              vy[offset] != std::numeric_limits<double>::max()) {
            sumx += vx[offset];
            sumy += vy[offset];
            cnt += 1;
          }
        }
      }
      if (cnt > 0) {
        vx[ind] = sumx / cnt;
        vy[ind] = sumy / cnt;
        holes[i] = -1;
      }
    }
    holes.erase(std::remove_if(holes.begin(), holes.end(), [](const int& x) {
        return x < 0;
    }), holes.end());
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
    for (int row=0; row!=nRows; ++row) {
      for (int col=0; col != nCols; ++col) {
        ux = *pvxForward;
        uy = *pvyForward;
        xp = col + t * ux;
        yp = row + t * uy;
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
    for (int row=0; row!=nRows; ++row) {
      for (int col=0; col != nCols; ++col) {
        ux = *pvxBackward;
        uy = *pvyBackward;
        xp = col + (1-t) * ux;
        yp = row + (1-t) * uy;
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
    std::vector<int> holes;
    findHoles(holes, vx, vy);
    fillHoles(holes, vx, vy);
    fillHoles(holes, vx, vy);
    killMaxLimits(holes, vx, vy);
}


void FlowInterpolation::colorTransfer(DImage& dest, const DImage& interp, const DImage& Im1, const DImage& Im2, const DImage& forward, const DImage& backward, double t) {

    const double *pInterp = interp.data();
    double *pDest = dest.data();
    const int nRows = Im1.height();
    const int nCols = Im1.width();
    double ux, uy, xp1, yp1, xp2, yp2;
    int xi1, yi1, xi2, yi2;

    if ( t > 0.5 ){
      dest.copyData(Im2);
    } else {
      dest.copyData(Im1);
    }
    dest.copyData(Im1);
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
