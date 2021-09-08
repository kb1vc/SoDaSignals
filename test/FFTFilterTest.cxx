#include "Filter.hxx"
#include <iostream>
#include "FFT.hxx"
#include <fstream>
#include <algorithm>
#include <cmath>
#include <complex>

std::ofstream debug("debug.out");
double magAtFreq1(double freq, std::vector<std::complex<double>> & ovec, SoDa::FFT & fft) {
  // take the FFT and shift
  std::vector<std::complex<double>> sv(ovec.size()); 
  fft.fft(sv, ovec);
  fft.shift(sv, sv);
  // find the mag at the index
  int ovsize = ovec.size();
  int idx = (ovsize >> 1) + (int) ((double (ovsize)) * freq);
  return std::abs(sv[idx]);

}

double magAtFreq(double freq, std::vector<std::complex<double>> & ovec, SoDa::FFT & fft) {
  // take the FFT and shift
  std::vector<std::complex<double>> sv(ovec.size());
  sv.resize(ovec.size());
  fft.fft(sv, ovec);
  fft.shift(sv, sv);
  // find the mag at the index
  int ovsize = ovec.size();
  int idx = (ovsize >> 1) + (int) ((double (ovsize)) * freq);


  return std::abs(sv[idx]);  
}

double magAtFreq2(double freq, std::vector<std::complex<double>> & ovec, SoDa::FFT & fft) {
  // take the FFT and shift
  std::vector<std::complex<double>> sv(ovec.size()); 
  fft.fft(sv, ovec);
  fft.shift(sv, sv);
  // find the mag at the index
  int ovsize = ovec.size();
  int idx = (ovsize >> 1) + (int) ((double (ovsize)) * freq);

  double ourpower = 0.0;
  for(int i = -2; i < 3; i++) {
    auto & v = sv[idx+i];
    
    ourpower = v.real() * v.real() + v.imag() * v.imag();
  }
  ourpower = ourpower / 5.0;

  // find the median power in a bin in the stop band 
  std::vector<double> spow(ovsize);
  for(int i = 0; i < ovsize; i++) {
    auto & v = sv[i]; 
    spow[i] = v.real() * v.real() + v.imag() * v.imag();
  }

  std::nth_element(spow.begin(), spow.begin() + (ovsize >> 2), spow.end());  
  double medpow = spow[ovsize >> 2];

  return ourpower / medpow;

}





int main() {
  int sample_rate = 1000; 
  int num_taps = 128;
  double f_sample_rate = ((double) sample_rate);
  int buf_len = 1000; 
  // a bandpass filter from 
  SoDa::Filter<double> filt_LP(SoDa::FilterType::BP, num_taps, 
		       f_sample_rate, 
		       0.0, 0.15 * f_sample_rate,
		       0.0025 * f_sample_rate,
			       50, buf_len);
  
  SoDa::Filter<double> filt_HP(SoDa::FilterType::BP, num_taps, 
		       f_sample_rate,
		       0.1 * f_sample_rate, 0.5 * f_sample_rate,
		       0.0025 * f_sample_rate,
			       50, buf_len);
  
  SoDa::Filter<double> filt_BP(SoDa::FilterType::BP, num_taps,
		       f_sample_rate,
		       -0.04 * f_sample_rate, 0.04 * f_sample_rate,
		       0.01 * f_sample_rate,
			       50, buf_len);
  
  SoDa::Filter<double> filt_BS(SoDa::FilterType::BS, num_taps, 
		       f_sample_rate,
		       -0.3 * f_sample_rate, -0.25 * f_sample_rate,
		       0.0025 * f_sample_rate,
			       50, buf_len);
  
  SoDa::FFT fft(buf_len);
  
  int buf_mult = 20; 
  for(double freq = -0.5 * f_sample_rate; 
      freq < 0.5 * f_sample_rate; 
      freq += 0.003125 * f_sample_rate) {
    // build a vector
    std::vector<std::complex<double>> test(buf_len);
    std::vector<std::complex<double>> out_LP(buf_len);
    std::vector<std::complex<double>> out_HP(buf_len);
    std::vector<std::complex<double>> out_BP(buf_len);
    std::vector<std::complex<double>> out_BS(buf_len);    
    double anginc = M_PI * 2.0 * freq / f_sample_rate; 

    for(int l = 0; l < buf_mult; l++) {
      double ang_start = anginc * ((double) (l * buf_mult));
      for(int i = 0; i < test.size(); i++) {
	double ang = ang_start + anginc * ((double) i); 
	test[i] = std::complex<double>(cos(ang), sin(ang)); 
      }
      
      // now apply the filter. 
      filt_LP.applyCont(out_LP, test);
      filt_HP.applyCont(out_HP, test);
      filt_BP.applyCont(out_BP, test);
      filt_BS.applyCont(out_BS, test);
    }

#if 1
    // now measure the magnitude of the last element of the output
    double mag_LP = std::abs(out_LP[out_LP.size() / 2]);
    double mag_HP = std::abs(out_HP[out_HP.size() / 2]);
    double mag_BP = std::abs(out_BP[out_BP.size() / 2]);
    double mag_BS = std::abs(out_BS[out_BS.size() / 2]);        
#else
    double mag_LP = magAtFreq(freq / f_sample_rate, out_LP, fft);
    double mag_HP = magAtFreq(freq / f_sample_rate, out_HP, fft);
    double mag_BP = magAtFreq(freq / f_sample_rate, out_BP, fft);
    double mag_BS = magAtFreq(freq / f_sample_rate, out_BS, fft);
#endif
#if 1    
    std::cout << freq << " "
	      << mag_HP << " "
	      << mag_LP << " "
	      << mag_BP << " "
	      << mag_BS << "\n";
#endif    
  }
}

