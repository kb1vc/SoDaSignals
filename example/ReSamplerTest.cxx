#include <SoDa/ReSampler.hxx>
#include <iostream>
#include <fstream>
#include <cmath>
#include <SoDa/Format.hxx>

int main() {
  int sample_rate = 1000; 
  // a bandpass filter from
  try {
    SoDa::ReSampler<float> resamp(500,
				  4,
				  5, 
				  0.05);

    std::cerr << "About to build test vectors\n";
    // build an input vector of 250 elements.
    // We're going to resample it to 200.
    std::vector<std::complex<float>> in(500), out(400);

    // put something in at 0.2 (f_nyquist / 2)

    float ang = 0.0;
    float ang_inc = 2 * M_PI / 23.0;

    std::ofstream inf("rs_in.dat");
    std::ofstream outf("rs_out.dat");

    std::cerr << "Building FFT\n";
    SoDa::FFT inxform(in.size());    
    std::vector<std::complex<float>> in_fft(in.size());
    SoDa::FFT outxform(out.size());
    std::vector<std::complex<float>> out_fft(out.size());

    int ki = 0;
    int ko = 0;

    std::cerr << "Creating sine\n";
    for(int tr = 0; tr < 5; tr++) {
      for(int i = 0; i < in.size(); i++) {
	in[i] = std::complex<float>(cos(ang), sin(-ang));
	std::cout << ang << " " << in[i].real() << " " << in[i].imag() << "\n";
	ang = ang + ang_inc;
	ang = (ang > M_PI) ? (ang - 2.0 * M_PI) : ang; 
      }

      std::cerr << "Resampling\n";

      std::cerr << " The input vector is getting fiddled with in some odd way. \n";
      resamp.apply(out, in);

      std::cerr << "Doing in/out xforms\n";
      // transform it
      inxform.fft(in_fft, in); 
      // transform the result
      outxform.fft(out_fft, out);

      std::cerr << "Printing\n";
      for(int i = 0; i < in.size(); i++, ki++) {
	inf << SoDa::Format("%0 %1 %2 %3 %4\n")
	  .addI(ki)
	  .addF(in[i].real())
	  .addF(in[i].imag())
	  .addF(in_fft[i].real())
	  .addF(in_fft[i].imag());
      }
    

      for(int i = 0; i < out.size(); i++, ko++) {
	outf << SoDa::Format("%0 %1 %2 %3 %4\n")
	  .addI(ko)
	  .addF(out[i].real())
	  .addF(out[i].imag())
	  .addF(out_fft[i].real())
	  .addF(out_fft[i].imag());
      }
    }

    inf.close();
    outf.close();
  }
  catch (std::runtime_error & e) {
    std::cerr << e.what() << "\n";
  }
  catch (...) {
    std::cerr << "I have no idea what just happened\n";
  }
}
