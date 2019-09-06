#include "SplatMotions.h"
#

void SplatMotions::splatMotionsBidirect(DImage& vx, DImage& vy, const DImage& vxForward, const DImage& vyForward, const DImage& vxBackward, const DImage& vyBackward, const DImage& Im1, const DImage& Im2, float t){


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
