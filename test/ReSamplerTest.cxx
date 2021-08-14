#include <ReSampler.hxx>
#include <iostream>
#include <fstream>
#include <cmath>
#include <SoDa/Format.hxx>

#define TEST_DTYPE float

int main() {
  int sample_rate = 1000; 
  // a bandpass filter from
  try {
    SoDa::ReSampler<TEST_DTYPE> resamp(200,
				  2, //2, // 4,
				  1, // 5, 
				  0.05);

    std::cerr << "About to build test vectors\n";
    // build an input vector of 250 elements.
    // We're going to resample it to 200.
    std::vector<std::complex<TEST_DTYPE>> in(200), out(400); // out(400);

    // put something in at 0.2 (f_nyquist / 2)

    TEST_DTYPE ang = 0.0;
    TEST_DTYPE ang_inc = 2 * M_PI / 23.0;

    std::ofstream inf("rs_in.dat");
    std::ofstream outf("rs_out.dat");

    std::cerr << "Building FFT\n";
    SoDa::FFT inxform(in.size());    
    std::vector<std::complex<TEST_DTYPE>> in_fft(in.size());
    SoDa::FFT outxform(out.size());
    std::vector<std::complex<TEST_DTYPE>> out_fft(out.size());

    int ki = 0;
    int ko = 0;

    std::cerr << "Creating sine\n";
    int num_trials = 1;
    for(int tr = 0; tr < num_trials; tr++) {
      for(int i = 0; i < in.size(); i++) {
	in[i] = std::complex<TEST_DTYPE>(cos(ang), sin(ang));
	//std::cout << "Hmmm\n"; // ang << " " << in[i].real() << " " << in[i].imag() << "\n";
	ang = ang + ang_inc;
	ang = (ang > M_PI) ? (ang - 2.0 * M_PI) : ang; 
      }

      std::cerr << "Resampling\n";

      std::cerr << " The input vector is getting fiddled with in some odd way. \n";
      resamp.applyCont(out, in);

      std::cerr << "Doing in/out xforms\n";
      // transform it
      inxform.fft(in_fft, in);
      //inxform.shift(in_fft, in_fft);
      // transform the result
      outxform.fft(out_fft, out);
      //outxform.shift(out_fft, out_fft);      

      std::cerr << "Printing\n";
      for(int i = 0; i < in.size(); i++, ki++) {
	auto v = in_fft[i];
	TEST_DTYPE m = v.real() * v.real() + v.imag() * v.imag();
	inf << SoDa::Format("%0 %1 %2 %3 %4 %5\n")
	  .addI(ki)
	  .addF(in[i].real())
	  .addF(in[i].imag())
	  .addF(in_fft[i].real())
	  .addF(in_fft[i].imag())
	  .addF(10.0 * log10(m));	  
      }
    

      for(int i = 0; i < out.size(); i++, ko++) {
	auto v = out_fft[i];
	TEST_DTYPE m = v.real() * v.real() + v.imag() * v.imag();
	
	outf << SoDa::Format("%0 %1 %2 %3 %4 %5\n")
	  .addI(ko)
	  .addF(out[i].real())
	  .addF(out[i].imag())
	  .addF(out_fft[i].real())
	  .addF(out_fft[i].imag())
	  .addF(10.0 * log10(m));
      }
    }

    std::cerr << "closing inf\n";
    inf.close();
    std::cerr << "closing outf\n";
    outf.close();
    std::cerr << "Releasing stuff\n";
  }
  catch (std::runtime_error & e) {
    std::cerr << "OOOPS!  ";
    std::cerr << e.what() << "\n";
  }
  catch (...) {
    std::cerr << "I have no idea what just happened\n";
  }
}
