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

TEST(Derivative, DerivativeDxAdvanceFilter) {
  int w = 3;
  int h = 2;
  int c = 3;
  FImage a(w, h, c);
  for (int i = 0; i != a.nelements(); ++i) {
    a[i] = std::rand() % RANDOM_MAX;
    EXPECT_GE(255, a[i]);
    EXPECT_LE(0, a[i]);
  }
  FImage adx;
  a.dx(adx, true); 

  float xFilter[5]={1,-8,0,8,-1};
  for(int i=0;i<5;i++)
      xFilter[i]/=12;
  FImage adx_tmp(w, h, c);

  float* pBuffer;
  float v;
  for (int i=0; i<h; i++) {
    for (int j=0; j<w; j++) {
      int offset = i * w * c;
      pBuffer = adx_tmp.data() + offset + j*c;
      for (int l=-2; l<=2; l++) {
        v = xFilter[l+2];
        int jj = __min(__max(j+l,0),w-1);
        for (int k=0; k<c; k++)
          pBuffer[k]+=a[offset+jj*c+k]*v;
      }
    }
  }

  for (int i = 0; i != a.nelements(); ++i) {
    EXPECT_NEAR(adx[i], adx_tmp[i], 1e-4);
  }
}

TEST(Derivative, DerivativeDy) {
  int w = 256;
  int h = 256;
  int c = 3;
  FImage a(w, h, c);
  for (int i = 0; i != a.nelements(); ++i) {
    a[i] = std::rand() % RANDOM_MAX;
    EXPECT_GE(255, a[i]);
    EXPECT_LE(0, a[i]);
  }
  FImage ady;
  a.dy(ady, false);

  FImage ady_tmp(w, h, c);
  for (int i = 0; i < h-1; i++) {
    for (int j = 0; j < w; j++) {
      int offset = i * w + j;
      for (int k = 0; k < c; k++) {
        ady_tmp[offset*c+k] = a[(offset+w)*c+k]-a[offset*c+k];
      }
    }
  }

  for (int i = 0; i != a.nelements(); ++i) {
    EXPECT_EQ(ady[i], ady_tmp[i]);
  }
}

TEST(Derivative, DerivativeDyAdvanceFilter) {
  int w = 3;
  int h = 2;
  int c = 3;
  FImage a(w, h, c);
  for (int i = 0; i != a.nelements(); ++i) {
    a[i] = std::rand() % RANDOM_MAX;
    EXPECT_GE(255, a[i]);
    EXPECT_LE(0, a[i]);
  }
  FImage ady;
  a.dy(ady, true); 

  float yFilter[5]={1,-8,0,8,-1};
  for(int i=0;i<5;i++)
      yFilter[i]/=12;
  FImage ady_tmp(w, h, c);

  float* pBuffer;
  float v;
  for (int i=0; i<h; i++) {
    for (int j=0; j<w; j++) {
      int offset = i * w * c;
      pBuffer = ady_tmp.data() + (i*w+j)*c;
      for (int l=-2; l<=2; l++) {
        v = yFilter[l+2];
        int ii = __min(__max(i+l,0),h-1);
        for (int k=0; k<c; k++)
          pBuffer[k]+=a[(ii*w+j)*c+k]*v;
      }
    }
  }

  for (int i = 0; i != a.nelements(); ++i) {
    EXPECT_NEAR(ady[i], ady_tmp[i], 1e-4);
  }
}

}
