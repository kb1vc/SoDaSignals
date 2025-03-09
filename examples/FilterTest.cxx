#include <SoDa/OSFilter.hxx>
#include <iostream>
#include <fstream>
#include <cmath>



int main() {
  int sample_rate = 1000; 
  // a bandpass filter from
  try {
    int sample_rate = 1000; 
    int num_taps = 33;
    int buflen = 400;    
    double f_sample_rate = ((double) sample_rate);
    // a bandpass filter from 
  SoDa::OSFilter filt_LP(0.0, // low cutoff
		       0.15 * f_sample_rate, // hi cutoff
		       0.0025 * f_sample_rate,  // skirt size
		       f_sample_rate, // sample rate
		       buflen); 


    std::vector<std::complex<float>> in(buflen), out(buflen);

    // put something in the passband

    float ang = 0.0;
    float ang_inc = 2 * M_PI * 0.005;

    int ko = 0;
    
    std::ofstream outf("ft_out.dat");

    for(int tr = 0; tr < 5; tr++) {
      for(int i = 0; i < in.size(); i++) {
	in[i] = std::complex<float>(cos(ang), sin(-ang));
	ang = ang + ang_inc;
	ang = (ang > M_PI) ? (ang - 2.0 * M_PI) : ang; 
      }
    

      // apply the filter
      filt_LP.apply(in, out);

      int lim = out.size() - 1;
      for(int i = 0; i < out.size(); i++, ko++) {
	float iv = (i == lim) ? (tr + 1.0) : 0;
	outf << " " << ko
	     << " " << in[i].real()
	     << " " << in[i].imag()
	     << " " << out[i].real()
	     << " " << out[i].imag()
	     << " " << iv
	     << "\n";
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
