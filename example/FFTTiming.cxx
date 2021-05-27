#include <SoDa/FFT.hxx>
#include <SoDa/Format.hxx>
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

int main(int argc, char * argv[]) {

  // trials, size 
  int trials = atof(argv[1]);
  int vsize = atof(argv[2]);

  std::vector<std::vector<std::complex<float>>> in_vecs, out_vecs; 

  std::vector<float> ai = { 0.1, 0.13242 };
  std::vector<float> ang = {0.0, 1.2}; 
  // create input buffers and output buffers
  for(int i = 0; i < trials; i++) {
    std::vector<std::complex<float>> iv(vsize);
    std::vector<std::complex<float>> ov(vsize);

    for(int j = 0; j < vsize; j++) {
      iv[j] = std::complex<float>(0.0, 0.0); 
      for(int aii = 0; aii < 2; aii++) {
	ang[aii] = bump_ang(ang[aii], ai[aii]);
	iv[j] += std::complex<float>(sin(ang[aii]), cos(ang[aii]));
      }
    }
    in_vecs.push_back(iv);
    out_vecs.push_back(ov);     
  }
  
  // now create the FFT object
  SoDa::FFT fft(vsize);

  // do the time trial
  auto start = std::chrono::high_resolution_clock::now();

  float sum = 0.0; 
  for(int t = 0; t < trials; t++) {
    fft.fft(out_vecs[t], in_vecs[t]);
    fft.ifft(in_vecs[t], out_vecs[t]); 
    sum += in_vecs[t][t].real();
  }

  auto end = std::chrono::high_resolution_clock::now();
  
  auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

  
  // how many points?
  double points = 2.0 * ((double) vsize) * ((double) trials);
  double elapsed = ns;

  //  std::cout << "ns " << ns << " elapsed " << elapsed << " sum = " << sum << "\n";
  // vsize, trials, elapsed, time-per-pt
  std::cout << SoDa::Format("%0 %1 %2 %3\n")
    .addI(vsize)
    .addI(trials)
    .addF(elapsed * 1e-9, 'e')
    .addF((elapsed / points) * 1e-9, 'e');
  
  
}
