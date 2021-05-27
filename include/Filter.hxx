#pragma once
#include <complex>
#include <vector>
#include "FFT.hxx"
#include "ChebyshevWindow.hxx"
#include <iostream>
#include <memory>

namespace SoDa {
  /**
   * \enum FilterType
   * 
   * Filter type -- band pass, or band stop. 
   */
  enum FilterType { BP, /*!< band pass */
		    BS  /*!< band stop */
  };
  
  template <typename T>
  class Filter {
  public:
    
    /**
     * @brief Constructor
     * 
     * All frequencies are in Hz, if sample_rate is in samples/second. 
     * 
     * 
     * 
     * @param typ Specify low pass, high pass, band pass, band stop filter
     * @param _num_taps The number of taps in the filter. (filter length)
     * @param sample_rate All frequencies are normalized to -0.5 sample_rate to 0.5 sample_rate
     * @param f1 The lower pass/stop band edge.  
     * @param f2 The upper pass/stop band edge
     * @param transition_width  Width of skirts (more or less)
     * @param stopband_atten Attenuation in dB for frequencies in the stop band, beyond the 
     * skirts. That's the magic of the Dolph-Chebyshev window. 
     *
     */
    Filter(FilterType typ, int _num_taps, float sample_rate, float f1, float f2, 
	   float transition_width, float stopband_atten)  {
      initFilter(typ, _num_taps, sample_rate, f1, f2, transition_width, stopband_atten); 
    }
    
  protected:
    Filter() {
      // do nothing constructor -- used for ReSampler, f'rinstance. 
    }
    void initFilter(FilterType typ, int _num_taps, float sample_rate, float f1, float f2, 
		    float transition_width, float stopband_atten) {

      num_taps = _num_taps; 

      // We *have* our standards, you know.
      if(stopband_atten < 25.0) stopband_atten = 25.0;
  
      // set all the relevant sizes
      h.resize(num_taps);

      H.resize(num_taps);
      // we haven't done an overlap and save op yet. 
      last_invec_size = 0;
      // we don't need an fft widget for the input stuff yet.
      dft_p = nullptr; 

      // make the prototype
      switch(typ) {
      case BP:
	makePrototype(H, sample_rate, f1, f2, transition_width);        
	break;
      case BS:
	makePrototype(H, sample_rate, f1 + transition_width, f2 - transition_width, transition_width);
	// now invert the filter
	for(auto & HE : H) {
	  HE = std::complex<T>(1.0, 0.0) - HE;
	}
	break;
      }

      // take the DFT
      SoDa::FFT fft(h.size());
      fft.ifft(h, H);

      // window it
      std::vector<T> cwin(num_taps);
      SoDa::ChebyshevWindow(cwin, num_taps, stopband_atten);

      // apply gain correction
      T gain_corr = 1.0 / ((T) num_taps); 
      for(int i = 0; i < num_taps; i++) {
	h[i] = h[i] * cwin[i] * gain_corr;
      }
    }

  public:
    /**
     * @brief Destructor
     * 
     * Free up FFT object, if we have one. 
     */
    ~Filter() {
    }
    
    /** 
     * @brief Apply the filter to a single isolated buffer
     * 
     * @param out The output of the filter. 
     * @param in The input to the filter. 
     */
    void apply(std::vector<std::complex<T>> & out, std::vector<std::complex<T>> & in) {
      if(H.size() != in.size()) makeImage(in.size());
  
      // take the FFT of the input
      dft_p->fft(out, in);
      // now multiply by the filter image
      // and the gain correction
      float gain_corr = 1.0 / ((T) in.size());
      for(int i = 0; i < in.size(); i++) {
	temp[i] = out[i] * H[i] * gain_corr; 
      }

      // now the inverse fft
      dft_p->ifft(out, temp); 
    }

    /** 
     * @brief Apply the filter to a buffer in a sequence of buffers. 
     * 
     * This method implements an overlap-and-save filter, suitable for a continous 
     * signal stream. 
     * 
     * @param out The output of the filter. 
     * @param in The input to the filter. 
     */
    void applyCont(std::vector<std::complex<T>> & out, std::vector<std::complex<T>> & in)
    {
      int i, j;
      debug_count++;
  
      if(last_invec_size != in.size()) {
	setupOverlap(in.size(), 1); 
      }

      // put the new samples at the end of the saved buffer.
      for(i = 0; i < in.size(); i++) {
	saved[i + save_length] = in[i];
      }
  
      // take the FFT of the saved buffer
      dft_p->fft(saved_dft, saved);

      // multiply by the image
      for(i = 0; i < saved_dft.size(); i++) {
	saved_dft[i] = H[i] * saved_dft[i];
      }
  
      // take the inverse fft
      dft_p->ifft(temp, saved_dft); 

      // copy the result to the output
      T scale = 1.0 / ((T) temp.size());
      for(i = 0; i < out.size(); i++) {
	out[i] = temp[i] * scale; 
      }

      // save the last section of the input
      for(i = 0, j = (in.size() - save_length); i < save_length; i++, j++) {
	saved[i] = in[j]; 
      }
    }

  protected:
    void setupOverlap(int invec_size, int required_multiple)
    {
      debug_count = 0; 
  
      // how big should the transform be?
      // Make it a power of two or a power of two times 3 or 5.
      int target_min = invec_size + num_taps;
      int target_max = target_min * 3 / 2; 

      std::cout << "target_min = " << target_min << "\ntarget_max = " << target_max << "\n";
      std::cout << "invec_size = " << invec_size << "\n";
      std::cout << "required multiple = " << required_multiple << "\n";

      // find a "good" size.
      // this is bone-brained stupid, but has the advantage that
      // I can explain it easily -- anything that's larger than
      // 16M entries isn't for us anyway.
      int good_size = 0;
      for(int i = 1; i < target_max; i = i * 2) {
	for(int j = 1; j < 8; j += 2) {
	  if((i * j * required_multiple) > target_min) {
	    int cand = i * j * required_multiple;
	    std::cout << "trying " << cand << "\n";
	    if(((cand < good_size) 
		|| (good_size == 0))
	       && (cand < target_max)) good_size = cand; 
	  }
	}
      }

      if(good_size == 0) {
	std::stringstream ss;
	ss << "SoDa::Filter:\n\tCould not find a good overlap-save length for a " 
	   <<  num_taps << " tap filter and " 
	   << invec_size << "%1 entry vector\n";
	throw std::runtime_error(ss.str());
      }
  
      std::cout << "Got size = " << good_size << "\n";
      
      if((good_size - invec_size) > (invec_size / 3)) {
	std::stringstream ss;
	ss << "SoDa::Filter:\n\tGot a really bad overlap-save  length of " 
	   << good_size 
	   << " for a " << num_taps  << " tap filter and " 
	   << invec_size << " entry vector\n";
	throw std::runtime_error(ss.str());
      }

  

      // clear the saved buffer
      for(int i = 0; i < saved.size(); i++) {
	saved[i] = std::complex<T>(0.0, 0.0);
      }
  
      saved_dft.resize(good_size);
  
      save_length = good_size - invec_size; 
  
      last_invec_size = invec_size;

      saved.resize(good_size);
  
      makeImage(good_size);
    }
    

  protected:
    ///@{
    int debug_count;
    
    int num_taps; 

    std::vector<std::complex<T>> h;
    std::vector<std::complex<T>> H;
    
    std::vector<std::complex<T>> saved;

    std::vector<std::complex<T>> saved_dft;
    
    std::vector<std::complex<T>> temp;
    
    int last_invec_size;

    std::shared_ptr<SoDa::FFT> dft_p;
    
    int save_length;
    ///@}
    
    /**
     * @name Protected Member Functions -- shake well before using
     */
    ///@{
    void makeImage(int image_size)  {
      H.clear();

      H.resize(image_size);

      std::cerr << "making image of size " << image_size << "\n";
      
      // create a new FFT object.
      dft_p = std::make_shared<SoDa::FFT>(image_size);

      std::vector<std::complex<T>> h_padded(image_size);

      int i, j; 
      // create the filter image.
      for(i = 0; i < image_size; i++) {
	h_padded[i] = std::complex<T>(0.0, 0.0);
      }

      for(i = 0; i < (num_taps + 1) / 2; i++, j++) {
	h_padded[i] = h[i]; 
      }
      for(i = 0; i < (num_taps / 2); i++) {
	h_padded[h_padded.size() - 1 - i] = h[h.size() - 1 - i];
      }

      temp.resize(image_size);
  
      // now transform it
      dft_p->fft(H, h_padded); 
    }
    
    void makePrototype(std::vector<std::complex<T>> & PH, float sample_rate, float f1, float f2, 
		       float transition_width) {
      // Assume we're building a band pass filter
      float lo_stop = f1 - transition_width;
      float lo_pass = f1;
      float hi_pass = f2;
      float hi_stop = f2 + transition_width;

      int num_taps = PH.size();
      float bucket_width = sample_rate / ((float) num_taps);

      // correct for lower band edge -- no transition band
      if(f1 < ((-0.5 * sample_rate) + bucket_width)) {
	lo_stop = -sample_rate; // now irrelevant
	lo_pass = -0.5 * sample_rate; 
      }
      // correct for upper band edge -- no transition band
      if(f2 > ((0.5 * sample_rate) - bucket_width)) {
	hi_stop = sample_rate;
	hi_pass = 0.5 * sample_rate; 
      }
  
      for(int i = 0; i < num_taps; i++) {
	int f_idx = i - (num_taps / 2);
	float freq = bucket_width * ((float) f_idx); 
	float v; 
	if((freq <= lo_stop) || (freq >= hi_stop)) {
	  PH[i] = 0.0;
	}
	else if((freq >= lo_pass) && (freq <= hi_pass)) {
	  PH[i] = 1.0;
	}
	else if((freq > lo_stop) && (freq < lo_pass)) {
	  v = (freq - lo_stop) / (lo_pass - lo_stop);
	  PH[i] = v; 
	}
	else if((freq > hi_pass) && (freq < hi_stop)) {
	  v = 1.0 - (freq - hi_pass) / (hi_stop - hi_pass);
	  PH[i] = v; 
	}
	else {
	  // I have no idea... 
	}
      }

      FFT::shift(PH, PH);
    }

    
    float normalize(float freq, float sample_rate) {
      float res = freq / sample_rate;

      res = M_PI * res;

      if(res > M_PI) res = res - 2.0 * M_PI;
      if(res < -M_PI) res = res + 2.0 * M_PI; 
      return res; 
    }
    ///@}
  };

  typedef Filter<float> FilterF;
  typedef Filter<double> FilterD;
}
