#include "FFT.hxx"
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <chrono>


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

void doTest(int s) {
  // create the FFT object
  SoDa::FFT fft(s); 

  // now create the vectors
  std::vector<std::complex<float>> iv(s);
  std::vector<std::complex<float>> ov(s);  
  // fill the input vector with something
  for(auto & v : iv) {
    v = std::complex<float>(1.0, 0.2);
  }

  // now do some timing runs -- 1000 iterations.
  int iters = 10;
  
  std::chrono::steady_clock::time_point t_start = std::chrono::steady_clock::now();
  for(int i = 0; i < iters; i++) {
    fft.fft(ov, iv); 
  }
  std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();  

  double ti = double(std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count());
  double ts = double(s * iters); 
  
  double cost_per_point = ti / ts ;
  double cost_per_vec =  ti / double(iters);
  
  std::cout << s << " " << ti << " " << cost_per_point << " " << cost_per_vec << "\n";
  std::cout.flush();
}

int main() {
  int lim = 256 * 1024; 
  int a2, a3, a5, a7, a11, a13;
  int i2, i3, i5, i7, i11, i13; 
  for(a2 = 512, i2 = 9; a2 < lim; a2 = a2 * 2, i2++) {
    for(a3 = a2, i3 = 0; a3 < lim; a3 = a3 * 3, i3++) {
      for(a5 = a3, i5 = 0; a5 < lim; a5 = a5 * 5, i5++) {
	for(a7 = a5, i7 = 0; a7 < lim; a7 = a7 * 7, i7++) {
	  for(a11 = a7, i11 = 0; (a11 < lim) && (i11 < 2); a11 = a11 * 11, i11++) {
	    for(a13 = a11, i13 = 0; (a13 < lim) && (i13 < 2); a13 = a13 * 13, i13++) {	    
	      char p2mark = (a13 == a2) ? 'P' : 'X';
	      std::cout << p2mark 
			<< " " << i2 << " " << i3 << " " 
			<< i5 << " " << i7 << " "
			<< i11 << " " << i13 << " ";
	      doTest(a13);
	    }
	  }
	}
      }
    }
  }
}
