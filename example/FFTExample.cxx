#include <iostream>
#include <math.h>
#include <vector>
#include <complex>
#include <SoDa/FFT.hxx>


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
  for(int i = 0; i < N; i++) {
    std::cout << in[i] << " " << out[i] / out_norm << " " << rev[i] / rev_norm << "\n";
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
