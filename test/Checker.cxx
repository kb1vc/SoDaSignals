#include "Checker.hxx"
#include <iostream>
#include <fstream>
#include <SoDa/Format.hxx>
#include "../include/NCO.hxx"

/*
  Copyright (c) 2022 Matthew H. Reilly (kb1vc)
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
  Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in
  the documentation and/or other materials provided with the
  distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
namespace SoDa {
  typedef std::vector<std::complex<float>> CVec;
  typedef std::vector<float> FVec;
  typedef std::vector<double> DVec; 

  float fixAngle(float v) {
    while(v > M_PI) v = v - 2.0 * M_PI; 
    while(v < -M_PI) v = v + 2.0 * M_PI; 
    return v; 
  }


  void dumpPhase(const std::string & fn, const DVec freqs, const CVec & a, const CVec & b) {
    std::ofstream os(fn);

    for(int i = 0; i < freqs.size(); i++) {
      os << SoDa::Format("%0 %1 %2 %3 %4\n")
			 .addF(freqs[i], 'e')
			 .addF(std::abs(a[i]))
			 .addF(std::arg(a[i]))
			 .addF(std::abs(b[i]))
			 .addF(std::arg(b[i]))
			 ;
    }
    os.close();
  }
  
  void dump2CVec(std::ostream & of, const CVec & a, const CVec & b) {
    
    for(int i = 0; i < std::min(a.size(), b.size()); i++) {
      of << SoDa::Format("%0 %1 %2 %3 %4\n")
	.addI(i)
	.addF(a[i].real(), 'e')
	.addF(a[i].imag(), 'e')
	.addF(b[i].real(), 'e')
	.addF(b[i].imag(), 'e');
    }
  }

  void dump2CVec(const std::string & fname, const CVec & a, const CVec & b) {
    std::ofstream of(fname); 
    dump2CVec(of, a, b); 
    of.close();
  }
  
  
  void dumpFVec(std::ostream & of, const FVec & a) {
    for(auto v : a) {
      of << SoDa::Format("%0\n").addF(v, 'e');
    }
  }

  void dumpFVec(const std::string & fname, const FVec & a) {
    std::ofstream of(fname); 
    dumpFVec(of, a); 
    of.close();
  }

  double Checker::phase(double freq) {
    double phi = (freq / high_sample_rate) * (double(filter_length + 1)) * M_PI;
    return -1.0 * fixAngle(phi);
  }

  Checker::Checker(double sample_rate, uint32_t filter_length, 
		   double ripple_limit_dB, double stopband_atten_dB, 
		   double permissible_phase_error,
		   uint32_t buffer_length,
		   uint32_t freq_steps) {

    initialize(sample_rate, filter_length,
	       ripple_limit_dB, stopband_atten_dB,
	       permissible_phase_error,
	       buffer_length,
	       buffer_length,
	       sample_rate,
	       freq_steps);
  }
  
  Checker::Checker(double sample_rate, uint32_t filter_length, 
		   double ripple_limit_dB, double stopband_atten_dB, 
		   double permissible_phase_error,
		   uint32_t input_buffer_length,
		   uint32_t output_buffer_length,
		   double input_sample_rate,
		   uint32_t freq_steps) {

    initialize(sample_rate, filter_length,
	       ripple_limit_dB, stopband_atten_dB,
	       permissible_phase_error,
	       input_buffer_length,
	       output_buffer_length,
	       input_sample_rate,
	       freq_steps);
  }

  void Checker::initialize(double _output_sample_rate, 
			   uint32_t _filter_length, 
			   double _ripple_limit_dB, 
			   double _stopband_atten_dB, 
			   double _permissible_phase_error,
			   uint32_t _input_buffer_length,
			   uint32_t _output_buffer_length,
			   double _input_sample_rate,
			   uint32_t _freq_steps) 
  {

    output_sample_rate = _output_sample_rate;
    filter_length = _filter_length; 
    ripple_limit_dB = _ripple_limit_dB;
    stopband_gain_dB = -_stopband_atten_dB;
    permissible_phase_error = _permissible_phase_error;
    input_buffer_length = _input_buffer_length;
    output_buffer_length = _output_buffer_length;
    input_sample_rate = _input_sample_rate;
    high_sample_rate = (output_sample_rate > input_sample_rate) ? output_sample_rate : input_sample_rate;
    
    uint32_t freq_steps = _freq_steps;

    ref_nco.setSampleRate(output_sample_rate);
    input_nco.setSampleRate(input_sample_rate);
    
    test_passed = true;

    // setup the oscillators
    first_output_oscillators.resize(freq_steps);
    second_output_oscillators.resize(freq_steps);
    first_input_oscillators.resize(freq_steps);
    second_input_oscillators.resize(freq_steps);    
    frequencies.resize(freq_steps);
    double lower_sample_rate = (input_sample_rate < output_sample_rate) ? input_sample_rate : output_sample_rate;
    double freq_step_size = lower_sample_rate / double(freq_steps);
    std::ofstream phase_out("phase_plot.dat");

    for(int i = 0; i < freq_steps; i++) {
      frequencies[i] = double(i) * freq_step_size - 0.5 * lower_sample_rate;

      ref_nco.setFreq(frequencies[i]);
      first_output_oscillators[i].resize(output_buffer_length);
      ref_nco.get(first_output_oscillators[i]);
      second_output_oscillators[i].resize(output_buffer_length);
      ref_nco.get(second_output_oscillators[i]);

      phase_out << frequencies[i] << " " << phase(frequencies[i]) << "\n";

      input_nco.setFreq(frequencies[i]);
      first_input_oscillators[i].resize(input_buffer_length);
      input_nco.get(first_input_oscillators[i]);
      second_input_oscillators[i].resize(input_buffer_length);
      input_nco.get(second_input_oscillators[i]);      
    }
    phase_out.close();
  }



  uint32_t Checker::getNumFreqSteps() {
    return first_output_oscillators.size();
  }

  double Checker::getFreq(uint32_t idx) { return frequencies[idx]; }

  void Checker::checkResponse(uint32_t freq_step,
			      std::function<CheckRegion(double)> freqRegion, 
			      std::function<void(std::vector<std::complex<float>> &,
						 std::vector<std::complex<float>>& )> filt) {
    
    // run the filter
    // twice
    std::vector<std::complex<float>> test_out(output_buffer_length);
    std::vector<float> freq_resp(second_output_oscillators.size());
    
    filt(first_input_oscillators[freq_step], test_out);
    // call it again
    filt(second_input_oscillators[freq_step], test_out);

    auto check_region = freqRegion(frequencies[freq_step]);
    
    // std::cerr << SoDa::Format("test freq %0  idx %1 reg %2\n")
    //   .addF(frequencies[freq_step], 'e')
    //   .addI(freq_step)
    //   .addI(int(check_region))
    //   ; 

    // find the autocorrelation for this step
    auto corr = SoDa::correlate(test_out, second_output_oscillators[freq_step]);

    target_phase_shift = phase(frequencies[freq_step]);
    // do the gains in dB power
    float gain = 20.0 * std::log10(std::abs(corr));
    float phase_shift = fixAngle(std::arg(corr));
    float phase_diff = fabs(phase_shift - target_phase_shift);
    auto orig_phase_diff = phase_diff;
    phase_diff = fabs(fixAngle(phase_diff));


    
    // now calculate the correlation for each frequency
    if(check_region == PASS_BAND) {
      // is the attenuation/gain less than threshold?
      if((gain < -ripple_limit_dB) && (gain > ripple_limit_dB)) {
	// we're out of the passband gain.
	std::cerr << SoDa::Format("Passband Gain out-of-range Frequency %0 gain %1 ripple limit %2\n")
	  .addF(frequencies[freq_step], 'e')
	  .addF(gain)
	  .addF(ripple_limit_dB)
	  ;
	if(test_passed) {	  
	  dump2CVec("ripple_fail.dat", test_out, second_output_oscillators[freq_step]);	    
	}
	test_passed = false;
      }

      // std::cerr << SoDa::Format("PSHIFT freq %0 target %1 shift %2 phase_diff %3 orig_diff %4 error %5\n")
      // 	.addF(frequencies[freq_step], 'e')
      // 	.addF(target_phase_shift)
      // 	.addF(phase_shift)
      // 	.addF(phase_diff)
      // 	.addF(orig_phase_diff)
      // 	.addF(abs(phase_shift - target_phase_shift))
      // 	;

      if(phase_diff > permissible_phase_error) {
	// the delay through the filter is wrong. 
	// we're the first failure.
	std::cerr << SoDa::Format("Passband shift out-of-range Frequency %0 shift %1 (%3) target %2\n")
	  .addF(frequencies[freq_step], 'e')
	  .addF(phase_diff)
	  .addF(target_phase_shift)
	  .addF(orig_phase_diff)
	  ;

	dump2CVec(SoDa::Format("phase_fail%0.dat").addF(frequencies[freq_step]).str(),
		    test_out, second_output_oscillators[freq_step]);


	test_passed = false;
      }
      else if((frequencies[freq_step] > 12.5e3) && (frequencies[freq_step] < 12.7e3)) {
	dump2CVec(SoDa::Format("phase_good%0.dat").addF(frequencies[freq_step]).str(),
		    test_out, second_output_oscillators[freq_step]);
	
      }
      
    }  // end check region PASS_BAND and on the test frequency
    else if(check_region == STOP_BAND) {
      if(gain > stopband_gain_dB) {
	std::cerr << SoDa::Format("Stopband gain exceeds limit, Frequency %0 gain %1 target %2\n")
	  .addF(frequencies[freq_step], 'e')
	  .addF(gain)
	  .addF(stopband_gain_dB)
	  ;
	if(test_passed) { dump2CVec("stopband_fail.dat", test_out, second_output_oscillators[freq_step]); }
	test_passed = false; 	
      }
    }
    else {
      // we're on the skirts.
      // anything goes
    }
    
    
    return; 
  }
}

