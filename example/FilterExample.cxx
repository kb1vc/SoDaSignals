#include <SoDa/Filter.hxx>
#include <iostream>

int main() {
  int sample_rate = 1000; 
  int num_taps = 128;
  double f_sample_rate = ((double) sample_rate);
  // a bandpass filter from 
  SoDa::Filter<double> filt_LP(SoDa::FilterType::BP, num_taps, 
		       f_sample_rate, 
		       0.0, 0.15 * f_sample_rate,
		       0.0025 * f_sample_rate,
		       50, sample_rate);
  
  SoDa::Filter<double> filt_HP(SoDa::FilterType::BP, num_taps, 
		       f_sample_rate,
		       0.1 * f_sample_rate, 0.5 * f_sample_rate,
		       0.0025 * f_sample_rate,
		       50, sample_rate);
  
  SoDa::Filter<double> filt_BP(SoDa::FilterType::BP, num_taps,
		       f_sample_rate,
		       -0.04 * f_sample_rate, 0.04 * f_sample_rate,
		       0.01 * f_sample_rate,
		       50, sample_rate);
  
  SoDa::Filter<double> filt_BS(SoDa::FilterType::BS, num_taps, 
		       f_sample_rate,
		       -0.3 * f_sample_rate, -0.25 * f_sample_rate,
		       0.0025 * f_sample_rate,
		       50, sample_rate);
  

  
  int buf_mult = 20; 
  for(double freq = -0.5 * f_sample_rate; 
      freq < 0.5 * f_sample_rate; 
      freq += 0.00125 * f_sample_rate) {
    // build a vector
    std::vector<std::complex<double>> test(sample_rate);
    std::vector<std::complex<double>> out_LP(sample_rate);
    std::vector<std::complex<double>> out_HP(sample_rate);
    std::vector<std::complex<double>> out_BP(sample_rate);
    std::vector<std::complex<double>> out_BS(sample_rate);    
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

    // now measure the magnitude of the last element of the output
    double mag_LP = abs(out_LP[out_LP.size() / 2]);
    double mag_HP = abs(out_HP[out_HP.size() / 2]);
    double mag_BP = abs(out_BP[out_BP.size() / 2]);
    double mag_BS = abs(out_BS[out_BS.size() / 2]);        
    std::cout << freq << " " 
	      << mag_LP << " "
	      << mag_HP << " "
	      << mag_BP << " "
	      << mag_BS << "\n";
  }
}
