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

namespace SoDa {

  /**
   * Here's the recipe
calculate the cutoff indexes in a discrete frequency domain view with M buckets
Hproto = 000111....0000...0
Hproto_r = fftshift(Hproto)
hproto = ifft(HProto_r)
// center the impulse
hproto_r = ifftshift(hproto)
hs = hproto_r * win

// now stuff the big image
h.resize(buffer_size)
h[0:hs.size...] = hs
H = fft(h)
   */
  void dumpCFVec(const std::string & fname, std::vector<std::complex<float>> & fv) {
    std::ofstream of(fname); 
    for(auto v : fv) {
      of << SoDa::Format("%0 %1\n").addF(v.real(), 'e').addF(v.imag(), 'e');
    }
    of.close();
  }

  void dumpFVec(const std::string & fname, std::vector<float> & fv) {
    std::ofstream of(fname); 
    for(auto v : fv) {
      of << SoDa::Format("%0 %1\n").addF(v, 'e');
    }
    of.close();
  }
  
  Filter::Filter(float low_cutoff, float high_cutoff, float skirt, 
		 float sample_rate, 
		 unsigned int image_size, 
		 float gain) {
    
    FilterSpec fspec(sample_rate, low_cutoff, high_cutoff, skirt);
    
    makeFilter(fspec, image_size, gain);
  }

  Filter::Filter(FilterSpec & filter_spec, 
		 unsigned int image_size,
		 float gain) {
    makeFilter(filter_spec, image_size, gain);
  }

  Filter::Filter(std::vector<std::complex<float>> & Hproto, 
		 unsigned int fft_size, float gain, 
		 WindowChoice window_choice) {
    sample_rate = 1.0; 
    makeFilter(Hproto, Hproto.size(), fft_size, gain, window_choice);
  }
  
  void Filter::makeFilter(std::vector<std::complex<float>> Hproto, 
			  unsigned int num_taps, 
			  unsigned int _image_size, 
			  float gain, 
			  WindowChoice window_choice) {
    image_size = _image_size;
    std::cerr << "#A";
    std::vector<std::complex<float>> hproto(num_taps);    

    std::cerr << "Hproto size " << Hproto.size() << " hproto size " << hproto.size() << "\n";
    // now we've got a frequency domain prototype.
    // first shift it to fit the FFT picture
    FFT pfft(num_taps);
    std::cerr << "#A1";
    pfft.shift(Hproto, Hproto);
    std::cerr << "#B";

    auto f0 = std::unique_ptr<FFT>(new FFT(image_size));    
    // transform it to the time domain
    pfft.ifft(Hproto, hproto);
    std::cerr << "#C";
    // shift that back to 0 in the middle
    pfft.ishift(hproto, hproto);
    std::cerr << "#D";
    // now apply a window
    if(window_choice != NOWINDOW) {
      std::vector<float> window(num_taps);
      //hammingWindow(window);
      hannWindow(window);    // gives the best adjacent frequency rejection
      //blackmanWindow(window);
      for(int i = 0; i < num_taps; i++) {
	hproto[i] = hproto[i] * window[i];
      }
    }

    // so now we have the time domain prototype.
    auto f2 = std::unique_ptr<FFT>(new FFT(image_size));
    std::cerr << "#E";
    // embed it in the impulse response of the appropriate length
    h.resize(image_size);
    auto f21 = std::unique_ptr<FFT>(new FFT(image_size));
    std::cerr << "image size " << image_size << " num_taps " << num_taps << "\n";
    std::cerr << "#E1";

    /// blows up after this
    for(int i = 0; i < num_taps; i++) {
      h[i] = hproto[i];
    }
    // blows up before this
    auto f22 = std::unique_ptr<FFT>(new FFT(image_size));
    std::cerr << "#E2";

    // zero the rest
    for(int i = num_taps; i < image_size; i++) {
      h[i] = std::complex<float>(0.0, 0.0);
    }
    auto f3 = std::unique_ptr<FFT>(new FFT(image_size));            
    std::cerr << "#F";
    // now make the frequency domain filter
    
    fft = std::unique_ptr<FFT>(new FFT(image_size));
    std::cerr << "#G";    
    H.resize(image_size);
    
    fft->fft(h, H);
    std::cerr << "#H";
    // and now we have it. But we need to rescale by image_size
    
    // find the max value and use this as 0dB
    float max = 0.0; 
    for(auto v : H) {
      auto vm = std::abs(v); 
      max = std::max(vm, max); 
    }
    std::cerr << "#I";
    float scale = 1.0 / max;
    
    for(auto & v : H) {
      v = v * scale * gain;
    }
    std::cerr << "#J";    
  }
    
  
  void Filter::makeFilter(FilterSpec & filter_spec, 
			  unsigned int _image_size, 
			  float gain) {
    sample_rate = filter_spec.getSampleRate();
    image_size = _image_size;
    auto num_taps = filter_spec.getTaps();    

    // first create a frequency domain ideal filter:
    std::vector<std::complex<float>> Hproto(num_taps);

    filter_spec.fillHproto(Hproto);

    makeFilter(Hproto, num_taps, image_size, gain);
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
    if((in_buf.size() != image_size) || (out_buf.size() != image_size)) {
      throw BadBufferSize("apply", in_buf.size(), out_buf.size(), image_size); 
    }
    if(temp_buf.size() != image_size) {
      temp_buf.resize(image_size); 
    }

    std::vector<std::complex<float>> * tbuf_ptr;
    if(in_out_mode.time_in) {
      // first do the transform
      fft->fft(in_buf, temp_buf);
      tbuf_ptr = & temp_buf; 
    }
    else {
      // we're already taking a transform as input. 
      tbuf_ptr = & in_buf; 
    }
    
    float scale = 1.0 / float(H.size());
    if(in_out_mode.time_out) {
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
    if((in_buf.size() != image_size) || (out_buf.size() != image_size)) {
      throw BadBufferSize("apply", in_buf.size(), out_buf.size(), image_size); 
    }

    if(temp_in_buf.size() != image_size) {
      temp_in_buf.resize(image_size); 
    }
    if(temp_out_buf.size() != image_size) {
      temp_out_buf.resize(image_size); 
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
  
  void Filter::dump(std::ostream & os) {
    for(int i = 0; i < H.size(); i++) {
      os << SoDa::Format("FILTER %0 %1 %2 %3 %4\n")
	.addI(i)
	.addF(H[i].real())
	.addF(H[i].imag())
	.addF(h[i].real())
	.addF(h[i].imag());
    }
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
}
