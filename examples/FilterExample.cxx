#include <SoDa/OSFilter.hxx>
#include <iostream>

int main() {
  int sample_rate = 1000; 
  int buffer_size = int(sample_rate);
  double f_sample_rate = ((double) sample_rate);
  // a bandpass filter from 
  SoDa::OSFilter filt_LP(0.0, // low cutoff
		       0.15 * f_sample_rate, // hi cutoff
		       0.0025 * f_sample_rate,  // skirt size
		       f_sample_rate, // sample rate
		       buffer_size); 
  
  SoDa::OSFilter filt_HP(0.15 * f_sample_rate, // low cutoff
		       0.5 * f_sample_rate, // hi cutoff
		       0.0025 * f_sample_rate,  // skirt size
		       f_sample_rate, // sample rate
		       buffer_size); 
  
  
  SoDa::OSFilter filt_BP(-0.04 * f_sample_rate, // low cutoff
		       0.04 * f_sample_rate, // hi cutoff
		       0.01 * f_sample_rate,  // skirt size
		       f_sample_rate, // sample rate
		       buffer_size); 

  
  int buf_mult = 20; 
  for(double freq = -0.5 * f_sample_rate; 
      freq < 0.5 * f_sample_rate; 
      freq += 0.00125 * f_sample_rate) {
    // build a vector
    std::vector<std::complex<float>> test(sample_rate);
    std::vector<std::complex<float>> out_LP(sample_rate);
    std::vector<std::complex<float>> out_HP(sample_rate);
    std::vector<std::complex<float>> out_BP(sample_rate);

    double anginc = M_PI * 2.0 * freq / f_sample_rate; 

    for(int l = 0; l < buf_mult; l++) {
      double ang_start = anginc * ((double) (l * buf_mult));
      for(int i = 0; i < test.size(); i++) {
	double ang = ang_start + anginc * ((double) i); 
	test[i] = std::complex<float>(cos(ang), sin(ang)); 
      }
      
      // now apply the filter. 
      filt_LP.apply(test, out_LP);
      filt_HP.apply(test, out_HP);
      filt_BP.apply(test, out_BP);
    }

    // now measure the magnitude of the last element of the output
    double mag_LP = abs(out_LP[out_LP.size() / 2]);
    double mag_HP = abs(out_HP[out_HP.size() / 2]);
    double mag_BP = abs(out_BP[out_BP.size() / 2]);
    std::cout << freq << " " 
	      << mag_LP << " "
	      << mag_HP << " "
	      << mag_BP << "\n";
  }
}
