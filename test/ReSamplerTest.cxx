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

#include "../include/ReSampler.hxx"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <SoDa/Options.hxx>
#include <SoDa/Format.hxx>

#include "../test/Checker.hxx"

typedef std::vector<std::complex<float>> CVec;
typedef std::vector<float> FVec; 

// the test strategy:
//
// Put various frequencies into the resampler across the input band.
// then ensure that the output magnitude is at the right level and,
// if in the pass band, at the appropriate frequency.
//
// This uses the checker from the filter test

SoDa::Checker::CheckRegion getRegion(double freq, double flo, double fhi, double skirt) {
  SoDa::Checker::CheckRegion reg; 
  if((freq > flo + 0.95 * skirt) && (freq < (fhi - 0.95 * skirt))) {
    reg = SoDa::Checker::PASS_BAND; 
    // std::cerr << SoDa::Format("Test freq %0 flo %1 fhi %2 skirt %3\n")
    //   .addF(freq, 'e').addF(flo, 'e').addF(fhi, 'e').addF(skirt, 'e');
  }
  else if((freq < (flo - 1.2 * skirt)) || (freq > (fhi + 1.2 * skirt))) {
    reg = SoDa::Checker::STOP_BAND;       
    // std::cerr << SoDa::Format("Test freq %0 flo %1 fhi %2 skirt %3\n")
    //   .addF(freq, 'e').addF(flo, 'e').addF(fhi, 'e').addF(skirt, 'e');
  }
  else {
    reg = SoDa::Checker::TRANSITION_BAND;             
  }

  return reg; 
}




int main(int argc, char * argv[]) {
  int sample_rate = 1000; 
  // a bandpass filter from

  double fs_in;
  double fs_out; 
  
  SoDa::Options cmd;
  cmd.add(&fs_in, "fsin", 'i', 48e3, "interpolate")
    .add(&fs_out, "fsout", 'd', 8e3, "decimate")
    .addInfo("Test rational resampler with sweeping input.\n");

  if(!cmd.parse(argc, argv)) {
    std::cout << "Bad command line\nFAILED\n";
    exit(-1);
  }

  bool is_ok = true;

  typedef std::pair<double, double> dpair; 
  std::list<dpair> test_freqs = {
    dpair(48e3, 8e3),
    dpair(625e3, 48e3),
    dpair(120e3, 48e3),
    dpair(1.2e6, 48e3),
    dpair(1.25e6, 48e3),
    dpair(2.5e6, 48e3)
  };

  
  if(cmd.isPresent("fsin") && cmd.isPresent("fsout")) {
    test_freqs.clear();
    test_freqs.push_back(dpair(fs_in, fs_out));
  }

  bool passed = true;
  
  for(auto fp : test_freqs) {
    SoDa::ReSampler resamp(fp.first, fp.second, 0.05);
    
    SoDa::Checker chk(fp.second,
		      resamp.getFilterLength(),
		      1.0,
		      50.0,
		      0.1,
		      resamp.getInputBufferSize(),
		      resamp.getOutputBufferSize(),
		      fp.first,
		      1024);
    std::cerr << SoDa::Format("Test %0 -> %1\n").addF(fp.first, 'e').addF(fp.second, 'e');

    double smaller_rate = (fp.first < fp.second) ? fp.first : fp.second;
    double f_hi = smaller_rate * 0.5;
    double f_lo = -f_hi;
    double skirt = 0.1 * smaller_rate; 
    
    
    for(uint32_t i = 0; i < chk.getNumFreqSteps(); i++) {
      chk.checkResponse(i, [f_lo, f_hi, skirt](double f) { return getRegion(f, f_lo, f_hi, skirt); },
			[& resamp](std::vector<std::complex<float>> & in,
				   std::vector<std::complex<float>> & out) { resamp.apply(in, out); });
      if(!chk.testPassed()) {
	passed = false; 
      }
    }
  }

  if(passed) {
    std::cout << "PASSED\n";
  }
  else {
    std::cout << "FAILED\n";
  }
  
}
