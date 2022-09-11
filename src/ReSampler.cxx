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

#include "ReSampler.hxx"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

namespace SoDa {
  typedef std::vector<std::complex<float>> CVec;  
  static void dumpCVec(const std::string & fn, CVec a) {
    std::ofstream of(fn);
    for(int i = 0; i < a.size(); i++) {
      of  << SoDa::Format("%0 %1 %2\n")
	.addI(i)
	.addF(a.at(i).real(), 'e')
	.addF(a.at(i).imag(), 'e');
    }
    of.close();
  }
  
  static void calcMag(const std::string & st, std::vector<std::complex<float>> & ve) {
    float max = 0.0;
    float sumsq = 0.0; 
    for(auto v : ve) {
      auto m = std::abs(v);
      max = std::max(m, max); 
      sumsq += m * m; 
    }
  }

  static uint32_t getGCD(uint32_t a, uint32_t b) {
    if(b == 0) return a;
    else if(a < b) {
      return getGCD(b, a);
    }
    else {
      return getGCD(b, a % b);
    }
  }


  ReSampler::ReSampler(float FS_in,
		       float FS_out,
		       float time_span_min) {
    uint32_t i_fs_in = ((uint32_t) FS_in);
    uint32_t i_fs_out = ((uint32_t) FS_out);
    auto gcd = getGCD(i_fs_in, i_fs_out);

    U = uint32_t(i_fs_out / gcd);
    D = uint32_t(i_fs_in / gcd);

    // this all works for D > U.  But for U > D we end up pitch shifting. 
    
    // How big is the low pass filter? 
    double corner_factor = 0.45; 
    double passband = std::min(FS_in, FS_out);
    double cutoff = passband * corner_factor; 

    double skirt_proportion = FS_in / (passband * (0.5 - corner_factor));

    double supression = 60.0;
    // use fred harris's estimate.
    // transition band is 10 % of FS_out
    uint32_t num_taps = int(0.5 + skirt_proportion * supression / 22.0); 
    if((num_taps % 2) == 0) num_taps++;
    if(num_taps < 121) num_taps = 121;
    // std::cerr << SoDa::Format("ReSampler: first guess at numtaps = %0\n").addI(num_taps);
    
    // now find the input buffer size -- make it long enough to span time_span_min
    uint32_t min_in_samples = uint32_t(FS_in * time_span_min);
    uint32_t min_out_samples = (U * min_in_samples) / D;

    // std::cerr << SoDa::Format("ReSampler:: min_in_samples = %0 min_out_samples = %1\n").addI(min_in_samples).addI(min_out_samples); 
    while((min_in_samples < 1000) || (min_out_samples < 1000)) {
      // need to goose min_in_samples until it is big enough
      min_in_samples += 1000;
      min_out_samples = (U * min_in_samples) / D;
      // std::cerr << SoDa::Format("ReSampler:: adj: min_in_samples = %0 min_out_samples = %1\n").addI(min_in_samples).addI(min_out_samples);       
    }
    
    scale_factor = float(D) / float(U);
    
    auto k = (min_in_samples + D - 1) / D;
    Lx = k * D;
    Ly = k * U;

    // std::cerr << SoDa::Format("ReSampler:: Lx = %0 Ly = %1 k = %2\n")
    //   .addI(Lx)
    //   .addI(Ly)
    //   .addI(k);
    // setup the save buffer -- it is at least as long as the filter, and must
    // be a multiple of D.
    int savek = (num_taps + D - 1) / D;
    save_count = savek * D;

    // it must be longer than the filter by at least one. 
    if(save_count < (num_taps + 1)) save_count = save_count + D;
    
    // remember our discard
    discard_count = save_count * U / D;

    // finally, the save window is one less than the number of taps
    num_taps = save_count + 1;

    // do we really want to do an LPF?
    lpf_p = std::unique_ptr<SoDa::Filter>(new SoDa::Filter(-cutoff, cutoff, 0.015 * cutoff, FS_in, 
							     num_taps, Lx));

    
    // create the input and output buffers
    x.resize(Lx);
    y.resize(Ly);
    X.resize(Lx);
    Y.resize(Ly);
    
    // zero the input buffer since we use the
    // end of it for the save buffer. 
    for(auto & s : x) {
      s = std::complex<float>(0.0,0.0);
    }
    // zero the output Y vector, as we may be upsampling
    for(auto & s : Y) {
      s = std::complex<float>(0.0,0.0);
    }

    // create the input and output FFTs.
    in_fft_p = std::unique_ptr<SoDa::FFT>(new SoDa::FFT(Lx));
    out_fft_p = std::unique_ptr<SoDa::FFT>(new SoDa::FFT(Ly));    
    
    // and that's it.
    // std::cerr << SoDa::Format("lpf taps %0  Lx %1 save_count %2 discard_count %3 inputBufferSize %4 outputBufferSize %5 U %6 D %7 fs_in %8 fs_out %9\n")
    //   .addI(num_taps)
    //   .addI(Lx)
    //   .addI(save_count)
    //   .addI(discard_count)
    //   .addI(getInputBufferSize())
    //   .addI(getOutputBufferSize())
    //   .addI(U)
    //   .addI(D)
    //   .addF(FS_in, 'e')
    //   .addF(FS_out, 'e');
  }


  ReSampler::~ReSampler() {
  }

  uint32_t ReSampler::getInputBufferSize() {
    return Lx - save_count;
  }

  
  uint32_t ReSampler::getOutputBufferSize() {
    return Ly - discard_count;
  }
  
  int apcount = 0; 
  uint32_t ReSampler::apply(std::vector<std::complex<float>> & in,
			    std::vector<std::complex<float>> & out) {
    if(in.size() != getInputBufferSize()) {
      throw BadBufferSize("Input", in.size(), getInputBufferSize());
    }

    if(out.size() != getOutputBufferSize()) {
      throw BadBufferSize("Output", out.size(), getOutputBufferSize());
    }

    // first do the overlap-and-save thing.
    for(int i = 0; i < save_count; i++) {
      x.at(i) = x.at(x.size() - save_count + i);
    }
    for(int i = save_count; i < Lx; i++) {
      x.at(i) = in.at(i - save_count); 
    }
    
    // now do the FFT
    in_fft_p->fft(x, X);

    // apply the filter -- we must do this
    // even if we are upsampling, as a raw copy of
    // an FFT segment looks like a brick wall otherwise. 
    lpf_p->apply(X, X, InOutMode(false,false));
      
    // now load the output Y vector
    if(Y.size() < X.size()) {
      auto y_half_count = ((Ly + 1)/ 2);
      for(int i = 0; i < y_half_count - 1; i++) {
	Y.at(i) = X.at(i);	
	Y.at(Ly - 1 - i) = X.at(Lx - 1 - i);	
      }
      Y.at(y_half_count) = X.at(y_half_count);
    }
    else {
      // we are up sampling. Y gets half of X in the bottom, half in the top.
      for(int i = 0; i < Lx; i++) {
	if(i < Lx/2) {
	  // we're on the DC and above side.
	  Y.at(i) = X.at(i);
	}
	else {
	  auto yi = (Ly - 1) - (Lx - 1) + i;
	  Y.at((Ly - 1) -(Lx - 1) + i) = X.at(i);
	}
      }
    }

    // dumpCVec("Y.dat", Y);    

    // do the inverse FFT
    out_fft_p->ifft(Y, y);
      
    // and copy to the output
    for(int i = 0; i < getOutputBufferSize(); i++) {
      out.at(i) = y.at(i + discard_count);
    }

    apcount++; 
    
    // and that's it!
    return 0;
  }

  uint32_t ReSampler::apply(float * in,
			    float * out) {
    return 0;
  }
}
