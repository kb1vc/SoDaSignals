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

void timeOperation(SoDa::FFT & fft, 
		   std::vector<std::complex<float>> & iv,
		   std::vector<std::complex<float>> & ov, 
		   int iters, 
		   const std::string & name) {

  std::chrono::steady_clock::time_point t_start = std::chrono::steady_clock::now();
  for(int i = 0; i < iters; i++) {
    fft.fft(ov, iv); 
  }
  std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();  


  double ti = double(std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count());
  double s = double(iv.size() * iters); 
  std::cout << name << " iterations " << iters << " on vector size " << iv.size() << " " 
	    <<  std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count() 
	    << "us\n"
	    <<  (ti / s) << " per point\n"    ;
}

void doSingle(std::vector<std::complex<float>> & iv,
	      std::vector<std::complex<float>> & ov,
	      int iters) {
  
  SoDa::FFT fft(iv.size()); 
  
  timeOperation(fft, iv, ov, iters, "Single");
}

void doInterleaved(std::vector<std::complex<float>> & iv,
		   std::vector<std::complex<float>> & ov,
		   int U, 
		   int iters) {
  SoDa::FFT fft(iv.size() / U, FFTW_MEASURE | FFTW_UNALIGNED, 
		U, U);

  timeOperation(fft, iv, ov, iters, "Interleaved");
}

int main(int argc, char * argv[]) {

  int bsize = atoi(argv[1]);
  int upsample = atoi(argv[2]); 
  int iters = atoi(argv[3]);
  int bufsize = bsize * upsample; 
  std::vector<std::complex<float>> iv(bufsize);
  std::vector<std::complex<float>> ov(bufsize);
  std::vector<std::complex<float>> siv(bsize);
  std::vector<std::complex<float>> sov(bsize);

  //  doSingle(siv, sov, iters * upsample);

  doSingle(iv, ov, iters);  

  //  doInterleaved(iv, ov, upsample, iters);   
  
}
