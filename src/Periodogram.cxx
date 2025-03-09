/*
  Copyright (c) 2022, 2025 Matthew H. Reilly (kb1vc)
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

#include "Periodogram.hxx"
#include "Filter.hxx"
#include <algorithm>

namespace SoDa {
  Periodogram::Periodogram(unsigned int segment_length, 
			   float _alpha,
			   Filter::WindowChoice window_choice) : 
    segment_length(segment_length) {

    setAlpha(_alpha);
    
    fft_scale = 1.0 / double(segment_length);
    
    acc_buffer.resize(segment_length);
    input_save_buffer.resize(segment_length);
    fft_in_buffer.resize(segment_length);
    fft_out_buffer.resize(segment_length);    

    input_save_buffer_valid_count = 0;
    accumulation_count = 0; 
    
    // build the FFT
    fft_p = std::unique_ptr<FFT>(new FFT(segment_length));

    // create the window
    window.resize(segment_length);
    switch (window_choice) {
    case Filter::HAMMING:
      Filter::hammingWindow(window);
      break; 
    case Filter::BLACKMAN:
      Filter::blackmanWindow(window);
      break; 
    case Filter::HANN:
    default:
      Filter::hannWindow(window);
      break; 
    }
    
    clear();
  }

  void Periodogram::setAlpha(const float _alpha) {
    alpha = abs(_alpha);     

    if(abs(alpha) < 1e-20) {
      alpha = 0.1; 
    }
    beta = 1.0 - alpha;
  }

  void Periodogram::accumulate(const std::vector<std::complex<float>> & in) {
    size_t curpos = 0; 

    while(curpos != in.size()) {
      // fill the end of the save_buffer with the input buffer
      size_t max_fill = segment_length - input_save_buffer_valid_count; 
      auto stuff_len = std::min(max_fill, in.size() - curpos); 

      std::copy(in.begin() + curpos, in.begin() + curpos + stuff_len, 
		input_save_buffer.begin() + input_save_buffer_valid_count);
      input_save_buffer_valid_count += stuff_len; 
      curpos += stuff_len; 

      // is the save_buffer full? -- if not, return; 
      if(input_save_buffer_valid_count != segment_length) return; 
      
      // else
      // apply the window
      for(int i = 0; i < segment_length; i++) {
	fft_in_buffer[i] = input_save_buffer[i] * window[i]; 
      }

      // do the FFT
      fft_p->fft(fft_in_buffer, fft_out_buffer); 

      // accumulate
      for(int i = 0; i < segment_length; i++) {
	acc_buffer[i] = alpha * std::abs(fft_out_buffer[i]) * fft_scale + beta * acc_buffer[i];
	//	acc_buffer[i] = fft_out_buffer[i] +  acc_buffer[i]; 	
      }
      accumulation_count++; 

      // move the last half of the buffer to the first half
      std::copy(input_save_buffer.begin() + (segment_length >> 1), 
		input_save_buffer.end(), 
		input_save_buffer.begin()); 
      
      input_save_buffer_valid_count = segment_length >> 1;
    }
  }

  void Periodogram::get(std::vector<float> & res) const {
    res.resize(acc_buffer.size());
    // fft shift the acc buffer into the result buffer
    int len = acc_buffer.size(); 
    int half_idx = len / 2;
    for(int i = 0; i < half_idx; i++) {
      res[i] = acc_buffer[half_idx + i];
      res[half_idx + i] = acc_buffer[i];
    }
    return; 
  }
    
  float Periodogram::getScaleFactor() {
    if(alpha == 1.0) {
      return 1.0 / float(accumulation_count); 
    }
    else {
      return 1.0; 
    }
  }
  
  void Periodogram::clear() {
    input_save_buffer_valid_count = 0; 
    accumulation_count = 0; 
    for(auto & v : acc_buffer) {
      v = 0.0; 
    }
  }
}
