#include "Checker.hxx"
#include <iostream>
#include <fstream>
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

  float fixAngle(float v) {
    while(v > M_PI) v = v - 2.0 * M_PI; 
    while(v < -M_PI) v = v + 2.0 * M_PI; 
    return v; 
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
    double phi = (freq / sample_freq) * (double(num_taps) + 1) * M_PI;

    phi = fmod(phi, 2.0 * M_PI);
    if(phi > M_PI) {
      phi = phi - 2.0 * M_PI; 
    }

    return phi;
  }
  
  Checker::Checker(double sample_freq, uint32_t filter_length, 
		   double ripple_limit_dB, double stopband_atten_dB, 
		   double permissible_phase_error,
		   uint32_t buffer_length,
		   uint32_t num_taps, 
		   uint32_t freq_steps) :
    sample_freq(sample_freq), filter_length(filter_length), 
    ripple_limit_dB(ripple_limit_dB),
    stopband_gain_dB(-stopband_atten_dB),
    permissible_phase_error(permissible_phase_error),
    buffer_length(buffer_length),
    num_taps(num_taps) 
  {
    
    ref_nco.setSampleRate(sample_freq);

    
    test_passed = true;

    // setup the oscillators
    first_oscillators.resize(freq_steps);
    second_oscillators.resize(freq_steps);    
    frequencies.resize(freq_steps);
    double freq_step_size = sample_freq / double(freq_steps);
    for(int i = 0; i < freq_steps; i++) {
      frequencies[i] = double(i) * freq_step_size - 0.5 * sample_freq;
      ref_nco.setFreq(frequencies[i]);
      first_oscillators[i].resize(buffer_length);
      ref_nco.get(first_oscillators[i]);
      second_oscillators[i].resize(buffer_length);
      ref_nco.get(second_oscillators[i]);
    }
  }

  uint32_t Checker::getNumFreqSteps() {
    return first_oscillators.size();
  }

  double Checker::getFreq(uint32_t idx) { return frequencies[idx]; }
  
  void Checker::checkResponse(
		       uint32_t freq_step,
		       std::function<CheckRegion(double)> freqRegion, 
		       std::function<void(std::vector<std::complex<float>> &,
					  std::vector<std::complex<float>>& )> filt) {

    // run the filter
    // twice
    std::vector<std::complex<float>> test_out(buffer_length);
    
    filt(first_oscillators[freq_step], test_out);
    // call it again
    filt(second_oscillators[freq_step], test_out);

    auto check_region = freqRegion(frequencies[freq_step]);
    
    std::cerr << SoDa::Format("test freq %0  idx %1 reg %2\n")
      .addF(frequencies[freq_step], 'e')
      .addI(freq_step)
      .addI(int(check_region))
      ; 


    
    // now calculate the correlation for each frequency
    for(int i = 0; i < second_oscillators.size(); i++) {
      // do the correlation for this frequency
      auto corr = SoDa::correlate(test_out, second_oscillators[i]);
      float gain = 10.0 * std::log10(std::abs(corr));
      float phase_shift = std::arg(corr);
      if((check_region == PASS_BAND) && (i == freq_step)) {
	//	std::cerr << "BANG  " << frequencies[i] << " " <<  phase_shift << "\n";
	// if we're in the passband, and at the same
	// frequency as the input, then we need to check
	// the ripple. 
      
	// is the attenuation/gain less than threshold?

	if((gain < -ripple_limit_dB) && (gain > ripple_limit_dB)) {
	  // we're out of the passband gain.
	  if(1 || test_passed) {
	    // we're the first failure.
	    std::cerr << SoDa::Format("Passband Gain out-of-range Frequency %0 gain %1 ripple limit %2 index %3 \n")
	      .addF(frequencies[i], 'e')
	      .addF(10 * std::log10(gain))
	      .addF(ripple_limit_dB)
	      .addI(i)
	      ;
	    dump2CVec("ripple_fail.dat", test_out, second_oscillators[i]);	    
	  }

	  test_passed = false;
	}

	target_phase_shift = phase(frequencies[i]);
	
	if(abs(phase_shift - target_phase_shift) > permissible_phase_error) {
	  // the delay through the filter is wrong. 
	  if(1 || test_passed) {
	    // we're the first failure.
	    std::cerr << SoDa::Format("Passband shift out-of-range Frequency %0 shift %1 target %2 index %3 \n")
	      .addF(frequencies[i], 'e')
	      .addF(phase_shift)
	      .addF(target_phase_shift)
	      .addI(i)
	      ;
	      
	    dump2CVec("phase_fail.dat", test_out, second_oscillators[i]);	  	    
	  }

	  test_passed = false;
	}
      }  // end check region PASS_BAND and on the test frequency
      else if((check_region == PASS_BAND) ||
	      (check_region == STOP_BAND) || 
	      ((check_region == TRANSITION_BAND) && (i != freq_step))) {
	// either we're in the passband but at a different test frequency
	// or in the stop band. 
	// or in the skirt and on a different frequency
	// in all cases we should be at the stopband level. 
	if(gain > stopband_gain_dB) {
	  if(1 || test_passed) {
	    std::cerr << SoDa::Format("Stopband gain exceeds passband limit, Frequency %0 gain %1 target %2 index = %3 \n")
	      .addF(frequencies[i], 'e')
	      .addF(gain)
	      .addF(stopband_gain_dB)
	      .addI(i)
	      ;
	    dump2CVec("stopband_fail.dat", test_out, second_oscillators[i]);	  	  	    
	  }

	  test_passed = false; 
	}
      }
      else {
	// we're on the skirts.
	// anything goes
      }
    }
    return; 
  }
}

