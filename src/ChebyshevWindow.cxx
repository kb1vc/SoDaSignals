#include "ChebyshevWindow.hxx"
#include "FFT.hxx"
#include <complex>
#include <iostream>

// Borrowed heavily from the scipy signal chebwin source code that was
// apparently borrowed from the octave source code.
//
// The Harris paper "On the Use of Win dows for Harmonic Analysis with
// the Discrete Fourier Transform" is the go to survey for all things
// window. (Proceedings of the IEEE Vol 66, No. 1 -- January 1978)
// However, the discussion of the Chebyshev window around page
// 193 in the Proceedings is somewhat incomplete.  In particular,
// the description of the Chebyshev function T_N is incorrect.
//
// Lyons ("Understanding Digital Signal Processing"  third edition) discusses
// the Chebyshev window in the FIR filter design chapter.  His description
// is also incomplete.  And in fact, it leaves out mention of the Chebyshev
// function all together.
//
// In the end, the source code for the scipy.signal.windows.chebwin was most
// helpful. https://github.com/scipy/scipy/blob/master/scipy/signal/windows/windows.py
// It shares a lot of its organization with the Octave code for the same function,
// and might have some common heritage.
//
// The chebwin function apparently wants to specify attenuation by amplitude (/20) vs.
// power.  That is the major change here.  
//

static float Cheb(int n, float x) {
  float ret; 
  if(abs(x) > 1.0) {
    float n1 = -1.0;
    if((n % 2) == 0) n1 = 1.0;
    if(x > 1.0) n1 = 1; 
    ret = n1 * cosh(n * acosh(abs(x)));
  }
  else {
    ret =  cos(n * acos(x));
  }
  return ret; 
}

void SoDa::ChebyshevWindow(std::vector<float> & res, size_t N, float atten) {
  res.resize(N);

  float nf = (float) N;
  float anginc = M_PI / nf;
  std::vector<float> Wk(N);
  float beta = cosh((1.0 / (nf - 1.0)) * acosh(pow(10.0, atten / 10.0)));

  std::vector<std::complex<float>> CW(N), cw(N);
  for(int i = 0; i < N; i++) {
    float ang = anginc * ((float) i);
    float x = beta * cos(ang);
    Wk[i] = Cheb(N - 1, x);
    Wk[i] = Wk[i];
    CW[i] = std::complex<float>(Wk[i] / Wk[0], 0.0);
    if((N % 2) == 0) CW[i] = CW[i] * std::exp(std::complex<float>(0.0, ang));
  }

  // now do the inverse fft to get us the actual window
  SoDa::FFT fft(N);
  fft.ifft(cw, CW);
  for(int i = 0; i < N; i++) {
    res[i] = cw[i].real() / cw[0].real();
  }
}
