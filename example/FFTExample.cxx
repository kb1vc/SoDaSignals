#include <SoDaFormat/Format.hxx>
#include <iostream>
#include <math.h>
#include <vector>
#include <complex>
#include "FFT.hxx"


int main() {
  unsigned int N = 16;
  SoDa::FFT fft(N);

  std::vector<std::complex<double>> in(N), out(N), rev(N);
  
  double angi = 2.0 * M_PI * 4.0 / ((double) N);
  double ang = 0.0;
  double cw[] = {1.0, 0.468, 0.014, 0.003, 0.002, 0.00115, 0.000614, 0.000266, 0,
		 0.000266,  0.000614, 0.00115, 0.002, 0.003, 0.014, 0.468};
  
  for(int i = 0; i < N; i++) {
    ang = ang + angi; 
    //    in[i] = std::complex<double>(cos(ang), 0.0); 
    in[i] = std::complex<double>(cw[i], 0.0);     
  }

		 
  fft.ifft(out, in);
  fft.fft(rev, out);
  
  double out_norm = out[0].real();
  double rev_norm = rev[0].real();
  std::cout << SoDa::Format("out_norm %0  rev_norm %1\n").addF(out_norm).addF(rev_norm);
  for(int i = 0; i < N; i++) {
    std::cout << SoDa::Format("%0 + j %1 ---> %2 + j %3 ---> %4 + j %5\n")
      .addF(in[i].real()).addF(in[i].imag())
      .addF(out[i].real() / out_norm).addF(out[i].imag() / out_norm)
      .addF(rev[i].real() / rev_norm).addF(rev[i].imag() / rev_norm);
  }


  // now check the shift function.
  for(int i = 6; i < 8; i++) {
    std::vector<int> iv(i), ov(i); 
    for(int j = 0; j < i; j++) {
      iv[j] = j; 
    }
    fft.shift(iv, iv); 
    for(int j = 0; j < i; j++) {
      std::cout << j << " " << iv[j] << "\n";
    }
  }
  
}
