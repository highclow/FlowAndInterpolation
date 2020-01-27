#include "gtest/gtest.h"

#include <cstdlib>
#include "Image.h"

namespace {

const int RANDOM_MAX = 256;

TEST(ImageLoader, DefaultLoader) {
  FImage a;
  EXPECT_EQ(true, a.IsEmpty());
  EXPECT_EQ(0, a.width());
  EXPECT_EQ(0, a.height());
  EXPECT_EQ(0, a.nchannels());
  EXPECT_EQ(0, a.npixels());
  EXPECT_EQ(0, a.nelements());
  EXPECT_EQ(true, a.IsFloat());
  EXPECT_EQ(false, a.isDerivativeImage());
}

TEST(ImageLoader, DimensionLoader) {
  int w = 6;
  int h = 5;
  int c = 3;
  FImage a(w, h, c);
  EXPECT_EQ(false, a.IsEmpty());
  EXPECT_EQ(w, a.width());
  EXPECT_EQ(h, a.height());
  EXPECT_EQ(c, a.nchannels());
  EXPECT_EQ(w*h, a.npixels());
  EXPECT_EQ(w*h*c, a.nelements());
  EXPECT_EQ(true, a.IsFloat());
  EXPECT_EQ(false, a.isDerivativeImage());
  for (int i = 0; i != a.nelements(); ++i) {
    EXPECT_EQ(0, a.data()[i]);
  }
}

TEST(ImageLoader, CopyLoader) {
  float v = 2.5;
  int w = 6;
  int h = 5;
  int c = 3;
  FImage a(v, w, h, c);
  FImage b(a);
  EXPECT_EQ(false, b.IsEmpty());
  EXPECT_EQ(w, b.width());
  EXPECT_EQ(h, b.height());
  EXPECT_EQ(c, b.nchannels());
  EXPECT_EQ(w*h, b.npixels());
  EXPECT_EQ(w*h*c, b.nelements());
  EXPECT_EQ(true, b.IsFloat());
  EXPECT_EQ(false, b.isDerivativeImage());
  for (int i = 0; i != b.nelements(); ++i) {
    EXPECT_EQ(2.5, b.data()[i]);
  }
}

TEST(ImageLoader, OperatorLoader) {
  float v = 2.5;
  int w = 6;
  int h = 5;
  int c = 3;
  FImage a(v, w, h, c);
  FImage b = a;
  EXPECT_EQ(false, b.IsEmpty());
  EXPECT_EQ(w, b.width());
  EXPECT_EQ(h, b.height());
  EXPECT_EQ(c, b.nchannels());
  EXPECT_EQ(w*h, b.npixels());
  EXPECT_EQ(w*h*c, b.nelements());
  EXPECT_EQ(true, b.IsFloat());
  EXPECT_EQ(false, b.isDerivativeImage());
  for (int i = 0; i != b.nelements(); ++i) {
    EXPECT_EQ(2.5, b.data()[i]);
  }
}

TEST(Derivative, DerivativeDx) {
  int w = 256;
  int h = 256;
  int c = 3;
  FImage a(w, h, c);
  for (int i = 0; i != a.nelements(); ++i) {
    a[i] = std::rand() % RANDOM_MAX;
    EXPECT_GE(255, a[i]);
    EXPECT_LE(0, a[i]);
  }
  FImage adx;
  a.dx(adx, false);

  FImage adx_tmp(w, h, c);
  for (int i = 0; i < h; i++) {
    for (int j = 0; j < w-1; j++) {
      int offset = i * w + j;
      for (int k = 0; k < c; k++) {
        adx_tmp[offset*c+k] = a[(offset+1)*c+k]-a[offset*c+k];
      }
    }
  }

  for (int i = 0; i != a.nelements(); ++i) {
    EXPECT_EQ(adx[i], adx_tmp[i]);
  }
}

//TEST(Derivative, DerivativeDxAdvanceFilter) {
//  int w = 256;
//  int h = 256;
//  int c = 3;
//  FImage a(w, h, c);
//  for (int i = 0; i != a.nelements(); ++i) {
//    a[i] = std::rand() % RANDOM_MAX;
//    EXPECT_GE(255, a[i]);
//    EXPECT_LE(0, a[i]);
//  }
//  FImage adx;
//  a.dx(adx, true); 
//
//  FImage adx_tmp(w, h, c);
//  for (int i = 0; i < h; i++) {
//    for (int j = 0; j < w-1; j++) {
//      int offset = i * w + j;
//      for (int k = 0; k < c; k++) {
//        adx_tmp[offset*c+k] = a[(offset+1)*c+k]-a[offset*c+k];
//      }
//    }
//  }
//
//  for (int i = 0; i != a.nelements(); ++i) {
//    EXPECT_EQ(adx[i], adx_tmp[i]);
//  }
//}

}
