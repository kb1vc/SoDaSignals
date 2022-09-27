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

#include "OSFilter.hxx"
#include <iostream>
#include <fstream>
#include "FFT.hxx"

namespace SoDa {

  
  OSFilter::OSFilter(float low_cutoff, float high_cutoff, float skirt, 
		 float sample_rate, unsigned int buffer_size) {

    // taps will be something like 2 * sample_rate / skirt  + odd;
    uint32_t gtaps = int((2 * sample_rate / skirt));
    // add 4 and make gtaps even (for a reason.)
    gtaps = (gtaps + 4) & ~0x1;
    // find a good FFT size
    uint32_t good_size = FFT::findGoodSize(buffer_size + gtaps);
    // taps are what's left over
    uint32_t taps = buffer_size - good_size + 1;

    FilterSpec fspec(sample_rate, gtaps);
    fspec
      .add(low_cutoff, 1.0)
      .add(high_cutoff + skirt * 0.5, 0.0);

    fspec.setTaps(taps);
    makeOSFilter(fspec, buffer_size); 
  }

  OSFilter::OSFilter(FilterSpec & filter_spec, 
		 unsigned int buffer_size) {
    makeOSFilter(filter_spec, buffer_size); 
  }

  OSFilter::OSFilter() {
    // don't do anything.  We'll set it up in a little while. 
  }
  
  void OSFilter::makeGenericFilter(std::vector<std::complex<float>> & H,
				   unsigned int num_taps, 
				   unsigned int _buffer_size, 
				   float gain,
				   WindowChoice window_choice
				   ) {
    buffer_size = _buffer_size; 

    unsigned int fft_size = buffer_size + num_taps - 1;
    std::cerr << "makeGenericFilter fft_size = " << fft_size << "\n";
    filter_p = std::unique_ptr<Filter>(new Filter(H, fft_size, gain, window_choice));
    // now size all the buffers. 
    x_augmented.resize(fft_size);
    y_augmented.resize(fft_size);
    X.resize(fft_size);
    Y.resize(fft_size);
    save_buf.resize(num_taps - 1);
    real_in.resize(buffer_size);
    real_in.resize(buffer_size);

    std::cerr << SoDa::Format("fft_size %0 num_taps %1 buffer_size %2\n")
      .addI(fft_size)
      .addI(num_taps)
      .addI(buffer_size);
  }
  
  void OSFilter::makeOSFilter(FilterSpec & filter_spec, 
			      unsigned int _buffer_size) {

    buffer_size = _buffer_size; 

    // taps will be something like 2 * sample_rate / skirt  + odd;
    uint32_t gtaps = filter_spec.estimateTaps(50, 500);
    
    // add 4 and make gtaps even (for a reason.)
    gtaps += 4;
    gtaps = (gtaps | 1) ^ 1;
    
    // find a good FFT size
    uint32_t good_size = FFT::findGoodSize(buffer_size + gtaps);
    // taps are what's left over
    uint32_t taps = good_size - buffer_size + 1;
    std::cerr << SoDa::Format("Buffer_size %0 good_size %1 taps %2 gtaps %3\n")
      .addI(buffer_size)
      .addI(good_size)
      .addU(taps)
      .addU(gtaps);
    
    filter_spec.setTaps(taps);
    
    filter_p = std::unique_ptr<Filter>(new Filter(filter_spec, good_size));
    // now size all the buffers. 
    x_augmented.resize(good_size);
    y_augmented.resize(good_size);
    X.resize(good_size);
    Y.resize(good_size);
    save_buf.resize(taps - 1);
    real_in.resize(buffer_size);
    real_in.resize(buffer_size);
  }


  unsigned int OSFilter::apply(std::vector<std::complex<float>> & in_buf, 
			       std::vector<std::complex<float>> & out_buf) {
    if((in_buf.size() != buffer_size) || (out_buf.size() != buffer_size)) {
      throw BadBufferSize("apply", in_buf.size(), out_buf.size(), buffer_size); 
    }
    
    // copy from the end of x_augmented to the start
    uint32_t end_start = x_augmented.size() - save_buf.size();
    for(int i = 0; i < save_buf.size(); i++) {
      x_augmented.at(i) = x_augmented.at(end_start + i); 
    }
    // now fill the rest
    for(int i = save_buf.size(); i < x_augmented.size(); i++) {
      int j = i - save_buf.size();
      x_augmented.at(i) = in_buf.at(j); 
    }
    // apply the filter.
    filter_p->apply(x_augmented, y_augmented); 

    // throw away the early samples
    for(int i = 0;  i < out_buf.size(); i++) {
      out_buf.at(i) = y_augmented.at(i + save_buf.size());
    }

    return out_buf.size(); 
  }

  unsigned int OSFilter::apply(std::vector<float> & in_buf, 
			       std::vector<float> & out_buf) {
    if((in_buf.size() != buffer_size) || (out_buf.size() != buffer_size)) {
      throw BadBufferSize("apply", in_buf.size(), out_buf.size(), buffer_size); 
    }

    for(int i = 0; i < in_buf.size(); i++) {
      real_in.at(i) = std::complex<float>(in_buf.at(i), 0.0);
    }
    
    apply(real_in, real_out);
    
    for(int i = 0; i < out_buf.size(); i++) {
      out_buf.at(i) = real_out.at(i).real();
    }
    return in_buf.size();    
  }
}
