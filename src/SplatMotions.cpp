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
        std::cout << ux << " " << uy << " " << xp << " " << yp << std::endl;
        splat(pixels, xp, yp);

        for (int k=0; k != 4; ++k) {
          int col = k << 1;
          int row = col + 1;
          int offset = nCols * i + j;
          cin.get();
          if (checkIndices(nRows, nCols, pixels[col], pixels[row])) {
          }
        }


        ++pvxForward;
        ++pvyForward;
        pIm1 += nChannels;
        pIm2 += nChannels;
      }
    }

    


//    splats = splat(splatty, xp, yp)
//    for s in splats:
//        if check_indices(splatty, s[0], s[1]):
//              old = splatty[s[1]][s[0]]
//              if np.isnan(old[0]) or np.isnan(old[1]):
//                splatty[s[1]][s[0]] = motion
//              else:
//                old_ptc = follow_intensity(frame0, frame1,
//                                           old, s[1], s[0], t)
//                new_ptc = follow_intensity(frame0, frame1,
//                                           motion, s[1], s[0], t)
//                if (new_ptc < old_ptc):
//                    splatty[s[1]][s[0]] = motion
}

void SplatMotions::splatMotionsBidirect(DImage& vx, DImage& vy, const DImage& vxForward, const DImage& vyForward, const DImage& vxBackward, const DImage& vyBackward, const DImage& Im1, const DImage& Im2, double t){

    const int nRows = vxForward.height();
    const int nCols = vxForward.width();
    std::cout << nRows << " " << nCols << std::endl;
    DImage pts;
    pts.allocate(nRows, nCols, 1);
    splatForward(vx, vy, pts, vxForward, vyForward, Im1, Im2, t);

//    h = frame0.shape[0]
//    w = frame0.shape[1]
//    splatty = np.zeros_like(forward)
//    # Flows start out undefined
//    splatty[:] = np.NAN
//    for row in tqdm(range(h), position=True, desc="splatting forward"):
//        for col in range(w):
//            splat_forward(forward, frame0, frame1, splatty,
//                          row, col, t)
//    for row in tqdm(range(h), position=True, desc="splatting backward"):
//        for col in range(w):
//            splat_backward(backward, frame0, frame1, splatty,
//                           row, col, t)
//
//    # Fill holes twice to get rid of big holes
//    fill_holes(splatty)
//    fill_holes(splatty)
//    kill_nans(splatty)
  
}
