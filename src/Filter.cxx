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

#include "Filter.hxx"
#include <iostream>
#include <fstream>
#include <SoDa/Format.hxx>

namespace SoDa {
  
  Filter::Filter(float low_cutoff, float high_cutoff, float skirt, 
		 float sample_rate,
		 unsigned int num_taps,
		 unsigned int buffer_size, 
		 float gain,
		 WindowChoice window_choice) {    
    
    FilterSpec fspec(sample_rate, low_cutoff, high_cutoff, skirt);

    fspec.setTaps(num_taps);
    
    makeFilter(fspec, buffer_size, gain, window_choice);
  }

  Filter::Filter(FilterSpec & filter_spec, 
		 unsigned int buffer_size,
		 float gain,
		 WindowChoice window_choice) {
    makeFilter(filter_spec, buffer_size, gain, window_choice);
  }

  Filter::Filter(std::vector<std::complex<float>> & Hproto, 
		 unsigned int buffer_size, 
		 float gain, 
		 WindowChoice window_choice) {
    sample_rate = 1.0; 
    makeFilter(Hproto, buffer_size, gain, window_choice);
  }
  
  void Filter::makeFilter(std::vector<std::complex<float>> Hproto, 
			  unsigned int _buffer_size, 
			  float gain, 
			  WindowChoice window_choice) {
    buffer_size = _buffer_size;

    auto num_taps = Hproto.size();
    
    std::vector<std::complex<float>> hproto(num_taps);    

    // now we've got a frequency domain prototype.
    // first shift it to fit the FFT picture
    FFT pfft(num_taps);

    pfft.shift(Hproto, Hproto);

    // transform it to the time domain
    pfft.ifft(Hproto, hproto);

    // shift that back to 0 in the middle
    pfft.ishift(hproto, hproto);

    // now apply a window
    std::vector<float> window(num_taps);    
    bool apply_window = true; 
    switch (window_choice) {
    case HAMMING:
      hammingWindow(window);
      break;
    case HANN:
      hannWindow(window);
      break;
    case BLACKMAN:
      blackmanWindow(window);
      break; 
    default:
      apply_window = false; 
    }

    if(apply_window) {
      for(int i = 0; i < num_taps; i++) {
	hproto[i] = hproto[i] * window[i];
      }
    }


    // so now we have the time domain prototype.
    // embed it in the impulse response of the appropriate length
    h.resize(buffer_size);
    auto f21 = std::unique_ptr<FFT>(new FFT(buffer_size));

    for(int i = 0; i < num_taps; i++) {
      h[i] = hproto[i];
    }

    // zero the rest
    for(int i = num_taps; i < buffer_size; i++) {
      h[i] = std::complex<float>(0.0, 0.0);
    }

    // now make the frequency domain filter
    
    fft = std::unique_ptr<FFT>(new FFT(buffer_size));
    H.resize(buffer_size);
    
    fft->fft(h, H);
    // and now we have it. But we need to rescale by buffer_size
    
    // find the max value and use this as 0dB
    float max = 0.0; 
    for(auto v : H) {
      auto vm = std::abs(v); 
      max = std::max(vm, max); 
    }
    float scale = 1.0 / max;
    
    for(auto & v : H) {
      v = v * scale * gain;
    }
  }
    
  
  void Filter::makeFilter(FilterSpec & filter_spec, 
			  unsigned int _buffer_size, 
			  float gain,
			  WindowChoice window) {
    sample_rate = filter_spec.getSampleRate();
    buffer_size = _buffer_size;
    auto num_taps = filter_spec.getTaps();    

    // first create a frequency domain ideal filter:
    std::vector<std::complex<float>> Hproto(num_taps);

    filter_spec.fillHproto(Hproto);

    makeFilter(Hproto, buffer_size, gain, window);
  }

  void Filter::hammingWindow(std::vector<float> & w) {
    int M = w.size();
    
    float anginc = 2 * M_PI / ((float) (M - 1));
    for(int i = 0; i < M; i++) {
      w[i] = 0.54 - 0.46 * cos(anginc * (float(i)));
    }
  }

  void Filter::hannWindow(std::vector<float> & w) {
    int M = w.size();
    
    float anginc = 2 * M_PI / ((float) (M - 1));
    for(int i = 0; i < M; i++) {
      w[i] = 0.5 - 0.5 * cos(anginc * (float(i)));
    }
  }

  void Filter::blackmanWindow(std::vector<float> & w) {
    int M = w.size();
    
    float anginc = 2 * M_PI / ((float) (M - 1));
    for(int i = 0; i < M; i++) {
      w[i] = 0.42 - 0.5 * cos(anginc * (float(i))) 
	+ 0.08 * cos(anginc * float(i) * 2.0);
    }
  }
  
  
  unsigned int Filter::apply(std::vector<std::complex<float>> & in_buf, 
			     std::vector<std::complex<float>> & out_buf, 
			     InOutMode in_out_mode) {
    if((in_buf.size() != buffer_size) || (out_buf.size() != buffer_size)) {
      throw BadBufferSize("apply", in_buf.size(), out_buf.size(), buffer_size); 
    }
    if(temp_buf.size() != buffer_size) {
      temp_buf.resize(buffer_size); 
    }

    std::vector<std::complex<float>> * tbuf_ptr;
    if(in_out_mode.xform_in) {
      // first do the transform
      fft->fft(in_buf, temp_buf);
      tbuf_ptr = & temp_buf; 
    }
    else {
      // we're already taking a transform as input. 
      tbuf_ptr = & in_buf; 
    }
    
    float scale = 1.0 / float(H.size());
    if(in_out_mode.xform_out) {
      // now multiply
      for(int i = 0; i < tbuf_ptr->size(); i++) {
	(*tbuf_ptr)[i] = (*tbuf_ptr)[i] * H[i] * scale;
      }
      // invert
      fft->ifft(*tbuf_ptr, out_buf); 
    }
    else {
      // they want frequency output
      // multiply directly into output buffer
      for(int i = 0; i < out_buf.size(); i++) {
	out_buf[i] = (*tbuf_ptr)[i] * H[i] * scale;
      }
    }
    return in_buf.size();
  }

  unsigned int Filter::apply(std::vector<float> & in_buf, 
			     std::vector<float> & out_buf, 
			     InOutMode in_out_mode) {
    if((in_buf.size() != buffer_size) || (out_buf.size() != buffer_size)) {
      throw BadBufferSize("apply", in_buf.size(), out_buf.size(), buffer_size); 
    }

    if(temp_in_buf.size() != buffer_size) {
      temp_in_buf.resize(buffer_size); 
    }
    if(temp_out_buf.size() != buffer_size) {
      temp_out_buf.resize(buffer_size); 
    }
    // first fill a complex vector
    for(int i = 0; i < in_buf.size(); i++) {
      temp_in_buf[i] = std::complex<float>(in_buf[i], 0.0); 
    }
    // now apply 
    apply(temp_in_buf, temp_out_buf);

    for(int i = 0; i < out_buf.size(); i++) {
      out_buf[i] = temp_out_buf[i].real();
    }
    return in_buf.size();    
  }
  
  
  std::pair<float, float> Filter::getFilterEdges() {
    // scan from the bottom and top to find the first
    // H sample over 0.5
    int hi, lo;
    std::vector<float> Himg(H.size());
    int half_size = H.size() / 2;
    for(int i = 0; i < H.size(); i++) {
      int j; 
      if(i < half_size) {
	j = half_size + i;
      }
      else {
	j = i - half_size; 
      }
      Himg[i] = std::abs(H[j]);
    }
    
    for(int i = 0; i < Himg.size(); i++) {
      if(std::abs(Himg[i]) > 0.5) {
	lo = i; 
	break; 
      }
    }
    for(int i = H.size() - 1; i >= 0; i--) {
      if(std::abs(Himg[i]) > 0.5) {
	hi = i; 
	break; 
      }
    }
    

    float hsize = float(H.size());
    float f_lo = float(lo);
    float f_hi = float(hi);
    float step = sample_rate / hsize;
    float bottom = -0.5 * sample_rate;
    auto ret = std::pair<float, float>(bottom + f_lo * step, bottom + f_hi * step);
    return ret; 
  }

  Filter::BadBufferSize::BadBufferSize(const std::string & st, 
				       unsigned int in, 
				       unsigned int out, 
				       unsigned int req) :
	std::runtime_error(SoDa::Format("Filter::%3 input and output buffer sizes (%0 and %1) must be equal to %2\n")
			   .addI(in)
			   .addI(out)
			   .addI(req)
			   .addS(st)
			   .str()) { }
  
}
