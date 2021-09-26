#include <ReSampler.hxx>
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <SoDa/Format.hxx>
#include <SoDa/Options.hxx>


#define TEST_DTYPE float

int main(int argc, char * argv[]) {
  int sample_rate = 1000; 
  // a bandpass filter from

  int interp;
  int decimate;
  
  SoDa::Options cmd;
  cmd.add(&interp, "interp", 'i', 1, "interpolate")
    .add(&decimate, "decim", 'd', 4, "decimate")
    .addInfo("Test rational resampler with noise input.\n");

  cmd.parse(argc, argv);
  
  try {
    int inbuflen = 7500;
    SoDa::ReSampler<TEST_DTYPE> resamp(inbuflen,
				  interp, //2, // 4,
				  decimate, // 5, 
				  0.05);

    std::cerr << "About to build test vectors\n";
    // build an input vector of 250 elements.
    // We're going to resample it to 200.
    std::vector<std::complex<TEST_DTYPE>> in(inbuflen);
    std::vector<std::complex<TEST_DTYPE>> out(in.size() * interp / decimate);
    
    std::ofstream outf("rs_out.dat");
    std::ofstream inf("rs_in.dat");

    // random generator
    std::random_device rd; // random device for seed
    std::mt19937 mt(rd()); // setup seed
    std::uniform_real_distribution<float> dist(-1.0, 1.0);
    
    std::cerr << "Building FFT\n";
    SoDa::FFT inxform(in.size());    
    std::vector<std::complex<TEST_DTYPE>> in_fft(in.size());
    SoDa::FFT outxform(out.size());
    std::vector<std::complex<TEST_DTYPE>> out_fft(out.size());

    std::cerr << SoDa::Format("In size = %0 out size = %1 interp = %2 decimate = %3\n")
      .addI(in.size())
      .addI(out.size())
      .addI(interp)
      .addI(decimate);
    
    int ki = 0;
    int ko = 0;

    std::cerr << "Creating sine\n";
    int num_trials = 3;
    for(int tr = 0; tr < num_trials; tr++) {
      for(int i = 0; i < in.size(); i++) {
	in[i] = std::complex<float>(dist(mt), dist(mt));
      }

      std::cerr << "Resampling\n";

      std::cerr << " The input vector is getting fiddled with in some odd way. \n";
      resamp.applyCont(out, in);

      std::cerr << "Doing in/out xforms\n";
      // transform it
      inxform.fft(in_fft, in);
      inxform.shift(in_fft, in_fft);
      // transform the result
      outxform.fft(out_fft, out);
      outxform.shift(out_fft, out_fft);      

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
