#include <ReSampler.hxx>
#include <iostream>
#include <fstream>
#include <cmath>
#include <SoDa/Format.hxx>

int main() {
  int sample_rate = 1000; 
  // a bandpass filter from
  try {
    int sample_rate = 48000; 
    int num_taps = 400;
    int buflen = 10000;    
    double f_sample_rate = ((double) sample_rate);
    // so frequencies can go from -24000 to 24000 Hz
    // a bandpass filter from 
    SoDa::Filter<float> filt_LP(SoDa::FilterType::BP, num_taps, 
				 f_sample_rate, 
				0.0, 80.0,
				20.0,
				buflen);


    std::vector<std::complex<float>> in(buflen), out(buflen);

    // put something in at 0.2 (f_nyquist / 2)

    float ang = 0.0;
    float f_test = 80.0;
    float ang_inc = 2 * M_PI * f_test / f_sample_rate;

    int ko = 0;
    
    std::ofstream outf("ft_out.dat");

    for(int tr = 0; tr < 5; tr++) {
      for(int i = 0; i < in.size(); i++) {
	in[i] = std::complex<float>(cos(ang), sin(ang));
	ang = ang + ang_inc;
	ang = (ang > M_PI) ? (ang - 2.0 * M_PI) : ang; 
      }
    

      // apply the filter
      filt_LP.applyCont(out, in);

      int lim = out.size() - 1;
      for(int i = 0; i < out.size(); i++, ko++) {
	float iv = (i == lim) ? (tr + 1.0) : 0;
	outf << SoDa::Format("%0 %1 %2 %3 %4 %5\n")
	  .addI(ko)
	  .addF(in[i].real())
	  .addF(in[i].imag())
	  .addF(out[i].real())
	  .addF(out[i].imag())
	  .addF(iv);
      }
    }

    outf.close();
  }
  catch (std::runtime_error & e) {
    std::cerr << e.what() << "\n";
  }
  catch (...) {
    std::cerr << "I have no idea what just happened\n";
  }
}
