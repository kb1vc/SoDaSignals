#include "FFT.hxx"
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <chrono>

float bump_ang(float a, float ai) {
  a = a + ai;
  if(a > M_PI) a = a - 2*M_PI;

  return a;
}

int main() {

  int bufsize = 256;
  std::vector<std::complex<float>> iv(bufsize);
  std::vector<std::complex<float>> fv(bufsize);  
  std::vector<std::complex<float>> ov(bufsize);

  SoDa::FFT fft(bufsize);
  
  float f_inc = 0.25; 
  for(float freq = 0.1; freq < 3.0; freq += f_inc) {
    float magerr = 0.0;
    float ang = 0.0;
    float ang_inc = freq;
    for(int i = 0; i < bufsize; i++) {
      iv[i] = std::complex<float>(cos(ang), sin(ang));
      ang = ang + ang_inc; 
    }

    fft.fft(fv, iv);

    // now invert it
    fft.ifft(ov, fv);

    // now plot it
    for(int i = 0; i < bufsize; i++) {
      ov[i] = ov[i] / ((float) bufsize); 
      std::cout << i << " " << iv[i].real() << " " << ov[i].real() << " " << (ov[i].real() - iv[i].real()) << "\n";
    }
  }
  
}
