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

#include "../include/FFT.hxx"
#include "../include/Hilbert.hxx"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <SoDa/Options.hxx>
#include <SoDa/Format.hxx>

typedef std::vector<std::complex<float>> CVec;
typedef std::vector<float> FVec; 

// the test strategy:
// 
// Pick a frequency and phase "randomly" or programmatically, and
// transform a real cos(2 pi f t + phi) using the hilbert transformer.
// Then compare the resulting output to ensure that the phase of the
// resulting time domain signal is in quadrature, and related (how?)
// to the input phase.
// 
// Also check the output signal amplitude and measure its RMS error
// with respect to a "locked" NCO. 
//

float fixAngle(float in) {
  float a = in;
  while(a > M_PI) a = a - 2.0 * M_PI;
  while(a < -M_PI) a = a + 2.0 * M_PI;
  return a;
}

void getStats(std::vector<float> & in, float & mean, float & variance) {
  mean = 0.0;
  float meansq;
  for(auto v : in) {
    mean += v;
    meansq += (v * v);
  }
  
  float N = float(in.size());

  mean = mean / N;
  meansq = meansq / N;

  variance = std::abs(meansq - mean * mean);
    
  return; 
}

struct Result {
  float freq; 
  float phase_diff_mean;
  float phase_diff_var;
  float mag_mean;
  float mag_var;
}; 

Result testHilbert(SoDa::Hilbert & xform, double freq, double phase, unsigned int blocks) {
  Result res;
  
  unsigned int vlen = xform.getSize();
  unsigned int num_taps = xform.getTaps();
  
  std::vector<float> test_in(vlen);
  std::vector<std::complex<float>> test_out(vlen);
  
  float angle = phase;
  float angle_incr = freq * M_PI;
  std::vector<float> angles(vlen);
  unsigned int acount = 0; 
  // do a number of passes, then check the last one
  for(int b = 0; b < blocks; b++) {
    for(int i = 0; i < vlen; i++) {
      test_in[i] = cos(angle);
      angles[i] = angle;
      angle += angle_incr; 
      angle = fixAngle(angle);
    }
    
    // now do the transform
    xform.apply(test_in, test_out);
  }
  
  // first calculate the phase and mag vectors
  std::vector<float> out_phase_diff(vlen);
  std::vector<float> out_mag(vlen);
  for(int i = 0; i < vlen; i++) {
    out_phase_diff[i] = fixAngle(std::arg(test_out[i]) - angles[i]);
    out_mag[i] = std::abs(test_out[i]);
  }

  res.freq = freq;
  getStats(out_phase_diff, res.phase_diff_mean, res.phase_diff_var);
  getStats(out_mag, res.mag_mean, res.mag_var);

  if((res.phase_diff_var > 0.5)  && freq > 0.3) {
    std::ofstream of("error.dat");
    for(int i = 0; i < vlen; i++) {
      of << SoDa::Format("%0 %1 %2\n")
	.addF(out_phase_diff[i])
	.addF(angles[i])
	.addF(std::arg(test_out[i]));
    }
    of.close();
    std::cerr << "Ooops  look at error.dat\n";
    exit(-1);
  }
  return res;
}

int main(int argc, char * argv[]) {
  double freq;
  double phase; 
  unsigned int num_blocks; 
  unsigned int image_size; 
  unsigned int taps;
  
  SoDa::Options cmd;
  cmd.add(&freq, "freq", 'f', 0.3, "test frequency value between 0 and 1.0")
    .add(&phase, "phase", 'p', 0.0, "test phase in range -pi to pi")
    .add(&image_size, "size", 's', 1024u, "size of input block")
    .add(&taps, "taps", 't', 0u, "number of taps in hilbert \"filter\"")
    .add(&num_blocks, "blocks", 'b', 4u, "number of blocks to process")
    .addInfo("Test Hilbert transform widget");

  if(!cmd.parse(argc, argv)) {
    std::cout << "Bad command line\nFAILED\n";
    exit(-1);
  }

  if(taps == 0) {
    taps = (image_size / 8) | 1;
  }
  bool is_ok = true;

  // create the transformer
  std::cerr << "In hilbert transform test.  Image size " << image_size << "\n";
  SoDa::Hilbert xform(taps, image_size);
  
  if(cmd.isPresent("freq")) {
    // we'll do just one test
    Result res = testHilbert(xform, freq, phase, num_blocks);
    std::cout << SoDa::Format("%0 %1 %2  %3  %4\n")
      .addF(res.freq, 'e')
      .addF(res.mag_mean, 'e')
      .addF(res.mag_var, 'e')            
      .addF(res.phase_diff_mean, 'e')
      .addF(res.phase_diff_var, 'e');
  }
  else {
    for(freq = 0.0; freq < 0.999; freq += 0.001327) {
      Result res = testHilbert(xform, freq, phase, num_blocks);
      std::cout << SoDa::Format("%0 %1 %2  %3  %4\n")      
	.addF(res.freq, 'e')
	.addF(res.mag_mean, 'e')
	.addF(res.mag_var, 'e')            
	.addF(res.phase_diff_mean, 'e')
	.addF(res.phase_diff_var, 'e');
    }
  }
}

