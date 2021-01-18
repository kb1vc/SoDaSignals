#include <SoDa/Filter.hxx>
#include <iostream>

int main() {
  int sample_rate = 1000; 
  int num_taps = 128;
  float f_sample_rate = ((float) sample_rate);
  SoDa::Filter filt_LP(SoDa::Filter::LO, num_taps, f_sample_rate, 0.15 * f_sample_rate, 0.01 * f_sample_rate, 0.0025 * f_sample_rate, 50);
  SoDa::Filter filt_HP(SoDa::Filter::HI, num_taps, f_sample_rate, 0.1 * f_sample_rate, 0.01 * f_sample_rate, 0.0025 * f_sample_rate, 50);  
  SoDa::Filter filt_BP(SoDa::Filter::BAND, num_taps, f_sample_rate, -0.1 * f_sample_rate, 0.18 * f_sample_rate, 0.0025 * f_sample_rate, 50);  

  
  int buf_mult = 4; 
  for(float freq = -0.5 * f_sample_rate; freq < 0.5 * f_sample_rate; freq += 0.00125 * f_sample_rate) {
    // build a vector
    std::vector<std::complex<float>> test(sample_rate);
    std::vector<std::complex<float>> out_LP(sample_rate);
    std::vector<std::complex<float>> out_HP(sample_rate);
    std::vector<std::complex<float>> out_BP(sample_rate);
    float anginc = M_PI * 2.0 * freq / f_sample_rate; 

    for(int l = 0; l < buf_mult; l++) {
      float ang_start = anginc * ((float) (l * buf_mult));
      for(int i = 0; i < test.size(); i++) {
	float ang = ang_start + anginc * ((float) i); 
	test[i] = std::complex<float>(cos(ang), sin(ang)); 
      }
      
      // now apply the filter. 
      filt_LP.applyCont(out_LP, test);
      filt_HP.applyCont(out_HP, test);
      filt_BP.applyCont(out_BP, test);
    }

    // now measure the magnitude of the last element of the output
    float mag_LP = abs(out_LP[out_LP.size() / 2]);
    float mag_HP = abs(out_HP[out_HP.size() / 2]);
    float mag_BP = abs(out_BP[out_BP.size() / 2]);    
    std::cout << freq << " " 
	      << mag_LP << " "
	      << mag_HP << " "
	      << mag_BP << "\n";
  }
}
