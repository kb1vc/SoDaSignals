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

#include "Hilbert.hxx"
#include <iostream>
#include <fstream>
#include "FFT.hxx"

namespace SoDa {

  Hilbert::Hilbert(unsigned int min_taps, unsigned int in_buffer_size) : 
    in_buffer_size(in_buffer_size) {
    
    std::cerr << "In Hilbert::Hilbert in_buffer_size = " << in_buffer_size << "\n";
    // create the input (complex) vector
    temp_in_buf.resize(in_buffer_size);
    temp_out_buf.resize(in_buffer_size);    

    // Let's find a good number of taps... ask the FFT thingie
    // we want at least min_taps taps
    auto good_size = FFT::findGoodSize(in_buffer_size + min_taps);
    auto taps = good_size - in_buffer_size + 1;

    
    std::cerr << "good_size " << good_size << "\n";
    // make the "filter".
    std::vector<std::complex<float>> H(taps);
    // the "negative" frequency side is zero
    // the reset is 2
    // except for the DC terms, which are 1/2
    int halfsize = (taps / 2) + 1;

    std::vector<std::complex<float>> Hp(taps); 
    for(int i = 0; i < halfsize; i++)  {
      Hp[i] = std::complex<float>(4.0, 0.0);
    }
    for(int i = halfsize; i < taps; i++) {
      Hp[i] = std::complex<float>(0.0,0.0);
    }
    // do the correction on the DC terms
    Hp[0] = Hp[0] * 0.5f;
    Hp[halfsize] = Hp[halfsize] * 0.5f;
    
    // we need to pre-shift this, as the generic filter expects
    // prototypes to be in FFT order
    FFT pfft(taps);
    pfft.shift(Hp, Hp);
    std::cerr << SoDa::Format("num taps %0 (%3) in_buffer_size %1 good_size %2\n")
      .addI(taps)
      .addI(in_buffer_size)
      .addI(good_size)
      .addI(Hp.size());
    
    makeGenericFilter(Hp, Hp.size(), in_buffer_size);

    debug_of.open("h_fft.dat");
  }
  
  unsigned int Hilbert::apply(std::vector<float> & in_buf, 
			      std::vector<std::complex<float>> & out_buf) {

    if((in_buf.size() != in_buffer_size) || (out_buf.size() != in_buffer_size)) {
      throw BadBufferSize("apply", in_buf.size(), out_buf.size(), in_buffer_size); 
    }

    // first create the input vector
    for(int i = 0; i < in_buffer_size; i++) {
      temp_in_buf[i] = std::complex<float>(2.0 * in_buf[i], 0.0);
    }
    
    // now apply this in the overlap-and-save manner
    return OSFilter::apply(temp_in_buf, out_buf); 
  }
}
