#include <ReSampler.hxx>
#include <iostream>
#include <fstream>
#include <cmath>
#include <memory>
#include <SoDa/Format.hxx>


// for a given resample ratio, sweep from -Fs/2 original to Fs/2 original
// and check for spikes in the frequency domain. 

#define TEST_DTYPE double

class FreqTest {
public:

  const int num_trials = 10;
  
  FreqTest(int interp, int decimate, TEST_DTYPE sample_rate) :
    interp(interp), decimate(decimate), sample_rate(sample_rate)
  {
    // pick an input buffer length that is a multiple of the decimation
    // rate
    in_buffer_length = 1024 * decimate;    
    in_buffer.resize(in_buffer_length);

    out_buffer_length = in_buffer_length * interp / decimate;
    out_buffer.resize(out_buffer_length); 
    
    in_stream.resize(in_buffer_length * num_trials);
    out_stream.resize(out_buffer_length * num_trials);
    out_err.resize(out_buffer_length * num_trials);    
    
    resamp_p = 
      std::unique_ptr<SoDa::ReSampler<TEST_DTYPE>>(
						   new SoDa::ReSampler<TEST_DTYPE>(in_buffer_length,
										   interp,
										   decimate,
										   0.05)); 
    
  }



  void doTest(TEST_DTYPE freq) {
    // first make the signal stream
    auto in_ph_incr = createStream(freq); 
    
    // skip the first buffer
    for(int trial = 0; trial < num_trials; trial++) {
      fillInBuffer(trial);
      resamp_p->applyCont(out_buffer, in_buffer);
      saveOutBuffer(trial); 
    }

    // check the output vector to see that the
    // phase advance from sample to sample is "correct". 

    auto out_ph_incr = in_ph_incr * ((TEST_DTYPE) decimate) / ((TEST_DTYPE) interp);
    // start from a good bit in on the output vector
    int start_idx = out_stream.size() / (2 * num_trials);
    TEST_DTYPE cur_angle = std::arg(out_stream[start_idx-1]); 
    max_err = 0.0;
    min_err = 1000.0;
    int max_idx; 
    TEST_DTYPE sum_sq_err = 0.0; 
    TEST_DTYPE new_angle; 
    for(int i = start_idx; i < out_stream.size(); i++) {
      new_angle = std::arg(out_stream[i]);      
      auto pa = phaseAdvance(cur_angle, new_angle); 
      auto p_error = abs(piCorr(pa - out_ph_incr));
      if(p_error > max_err) {
	max_err = p_error;
	max_idx = i; 
      }
      if(p_error < min_err) min_err = p_error;
      if(p_error > abs(out_ph_incr * 8.0)) {
	std::cerr << SoDa::Format("p_error %3 cur_angle %0 new_angle %1 pa %2 OPI %8  freq %4 i %5 %6/%7\n")
	  .addF(cur_angle)
	  .addF(new_angle)
	  .addF(pa)
	  .addF(p_error)
	  .addF(freq)
	  .addI(i)
	  .addI(interp)
	  .addI(decimate)
	  .addF(out_ph_incr);
      }
      sum_sq_err += p_error * p_error;
      out_err[i] = piCorr(pa - out_ph_incr);
      cur_angle = new_angle;
    }
    
    //    std::cerr << "Max idx = " << max_idx << " Max err = " << max_err << "\n";
    rms_error = sqrt(sum_sq_err / ((TEST_DTYPE) (out_stream.size() - start_idx)));
    
    rms_error = abs(rms_error / out_ph_incr);
    min_err = abs(min_err / out_ph_incr);
    max_err = abs(max_err / out_ph_incr);
  }

  TEST_DTYPE piCorr(TEST_DTYPE a) {
    while (a > M_PI) a = a - 2.0 * M_PI;
    while (a < -M_PI) a = a + 2.0 * M_PI;
    return a;
  }

  void getResults(TEST_DTYPE & r_e, TEST_DTYPE & min_e, TEST_DTYPE & max_e) {
    r_e = rms_error;
    min_e = min_err;
    max_e = max_err;
  }

  TEST_DTYPE phaseAdvance(TEST_DTYPE angle, TEST_DTYPE new_angle) {
    TEST_DTYPE ret = new_angle - angle;
    if(ret > M_PI) ret = ret - 2.0 * M_PI;
    if(ret < -M_PI) ret = ret + 2.0 * M_PI;
    return ret; 
  }

  void fillInBuffer(int tr) {
    int i = tr * in_buffer_length;
    
    for(auto & v : in_buffer) {
      v = in_stream[i++];
    }
  }

  void saveOutBuffer(int tr) {
    int i = tr * out_buffer_length;
    
    for(auto & v : out_buffer) {
      out_stream[i++] = v;
    }
  }
  
  
  TEST_DTYPE createStream(TEST_DTYPE freq) {
    TEST_DTYPE angle = 0.0;
    TEST_DTYPE ph_incr = M_PI * freq / sample_rate; 
    
    for(auto & v : in_stream) {
      v = std::complex<TEST_DTYPE>(cos(angle), sin(angle));
      angle += ph_incr; 
    }
    
    return ph_incr; 
  }

  void dump(std::string inname, std::string outname){
    std::ofstream ibufstr(inname);
    std::ofstream obufstr(outname);

    for(auto & v : in_stream) {
      ibufstr << v.real() << " " << v.imag() << "\n";
    }
    for(int i = 0; i < out_stream.size(); i++) {
      auto v = out_stream[i];
      auto e = out_err[i];
      obufstr << v.real() << " " << v.imag() << " "
	      << e.real() << " " << e.imag() << "\n";
    }
    
    ibufstr.close();
    obufstr.close();
    
  }
  
  int in_buffer_length;
  int out_buffer_length;
  std::vector<std::complex<TEST_DTYPE>> in_buffer, in_stream;
  std::vector<std::complex<TEST_DTYPE>> out_buffer, out_stream, out_err;  

  TEST_DTYPE rms_error, min_err, max_err;
  
  std::unique_ptr<SoDa::ReSampler<TEST_DTYPE>> resamp_p;
  int interp, decimate;
  TEST_DTYPE sample_rate; 
}; 

int main() {
  TEST_DTYPE sample_rate = 1000; 
  for(int interp = 1; interp < 9; interp++) {
    for(int decim = 1; decim < 9; decim++) {
      try {
	FreqTest ft(interp, decim, sample_rate);
	TEST_DTYPE freq_lim = 499.0 * ((TEST_DTYPE)interp)/((TEST_DTYPE)decim);
	if(freq_lim > 499.0) freq_lim = 499.0;
	for(TEST_DTYPE freq = - freq_lim; freq < freq_lim; freq += 0.9) {
	  ft.doTest(freq);
	  TEST_DTYPE r_e, min_e, max_e;
	  ft.getResults(r_e, min_e, max_e); 
	  std::cout << SoDa::Format("%4 %5 %0 %1 %2 %3\n")
	    .addF(freq)
	    .addF(r_e)
	    .addF(min_e)
	    .addF(max_e)
	    .addI(interp)
	    .addI(decim);
	}
      }
      catch (std::runtime_error & e) {
	std::cerr << SoDa::Format("Ooops -- runtime error %0 interp = %1 decimate = %2\n")
	  .addS(e.what())
	  .addI(interp)
	  .addI(decim);
      }
      catch (...) {
	std::cerr << SoDa::Format("I have no idea what just happened  interp = %0 decimate = %1\n")
	  .addI(interp)
	  .addI(decim);
      }
    }
  }
}
