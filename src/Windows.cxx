#include "Windows.hxx"
#include "FFT.hxx"
#include <complex>
#include <iostream>
#include <SoDa/Format.hxx>
#include <fstream>


/*
BSD 2-Clause License

Copyright (c) 2022 Matt Reilly - kb1vc
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/



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
  if(fabs(x) > 1.0) {
    float n1 = -1.0;
    if((n % 2) == 0) n1 = 1.0;
    if(x > 1.0) n1 = 1; 
    ret = n1 * cosh(n * acosh(fabs(x)));
  }
  else {
    ret =  cos(n * acos(x));
  }
  return ret; 
}

template <typename T>
void _ChebyshevWindow(std::vector<T> & res, size_t N, float atten) {
  res.resize(N);

  T nf = (T) N;
  T anginc = M_PI / nf;
  std::vector<T> Wk(N);
  T beta = cosh((1.0 / (nf - 1.0)) * acosh(pow(10.0, atten / 10.0)));

  std::vector<std::complex<T>> CW(N), cw(N);
  for(int i = 0; i < N; i++) {
    T ang = anginc * ((T) i);
    T x = beta * cos(ang);
    Wk[i] = Cheb(N - 1, x);
    Wk[i] = Wk[i];
    CW[i] = std::complex<T>(Wk[i] / Wk[0], 0.0);
    if((N % 2) == 0) CW[i] = CW[i] * std::exp(std::complex<T>(0.0, ang));
  }

  // now do the inverse fft to get us the actual window
  SoDa::FFT fft(N);
  fft.ifft(cw, CW);
  
  for(int i = 0; i < N; i++) {
    res[i] = cw[i].real() / cw[0].real();
  }
}


void SoDa::ChebyshevWindow(std::vector<float> & res, size_t N, float atten) {
  _ChebyshevWindow<float>(res, N, atten);
}

void SoDa::ChebyshevWindow(std::vector<double> & res, size_t N, float atten) {
  _ChebyshevWindow<double>(res, N, atten);
}

//
// This looks SOOOO much better than the dolph-chebyshev window.  I don't
// know what I was thinking back then.


template <typename T>
void _HammingHannWindow(std::vector<T> & ret, size_t N, T a0) {
  ret.resize(N);
  std::vector<T> res(N);
  
  float a1 = 1 - a0;

  T anginc = 2.0 * M_PI / (T(N-1));

  for(int n = 0; n < N; n++) {
    T ang = (T(n) * anginc);
    res[n] = a0 - a1 * cos(ang);
  }

  // do the fft shift
  for(int i = 0; i < N/2; i++) {
    ret[i] = res.at(N/2 + i);
    ret.at(N - i - 1) = res.at(N/2 - i);
  }

}


void SoDa::HammingWindow(std::vector<float> & res, size_t N) {
  float a0 = 25.0 / 46.0; 
  _HammingHannWindow<float>(res, N, a0);
}

void SoDa::HammingWindow(std::vector<double> & res, size_t N) {
  double a0 = 25.0 / 46.0;   
  _HammingHannWindow<double>(res, N, 0.54);
}

void SoDa::HannWindow(std::vector<float> & res, size_t N) {
  float a0 = 0.5;
  _HammingHannWindow<float>(res, N, a0);
}

void SoDa::HannWindow(std::vector<double> & res, size_t N) {
  double a0 = 0.5;
  _HammingHannWindow<double>(res, N, 0.54);
}
