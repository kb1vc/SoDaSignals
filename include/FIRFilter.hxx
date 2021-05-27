#pragma once
#include <complex>
#include <vector>
#include "ChebyshevWindow.hxx"
#include "FilterType.hxx"
#include "CircularBuffer.hxx"
#include "FFT.hxx"
#include <iostream>

namespace SoDa {
  /**
   * @brief Time domain implementation of a FIR filter
   */
  template <typename T>
  class FIRFilter {
  public:
    
    /**
     * @brief Constructor
     * 
     * All frequencies are in Hz, if sample_rate is in samples/second. 
     * 
     * 
     * 
     * @param typ Specify low pass, high pass, band pass, band stop filter
     * @param num_taps The number of taps in the filter. (filter length)
     * @param sample_rate All frequencies are normalized to -0.5 sample_rate to 0.5 sample_rate
     * @param f1 The lower pass/stop band edge.  
     * @param f2 The upper pass/stop band edge
     * @param transition_width  Width of skirts (more or less)
     * @param stopband_atten Attenuation in dB for frequencies in the stop band, beyond the 
     * skirts. That's the magic of the Dolph-Chebyshev window. 
     *
     */
    FIRFilter(FilterType typ, int num_taps, float sample_rate, float f1, float f2, 
	      float transition_width, float stopband_atten)  :
      num_taps(num_taps) {

      // create a circular buffer for the inputs.  This allows us to
      // process the input signal in chunks.  The circular buffer is the
      // "Z register" or delay line.
      Z_buffer.resize(num_taps, T(0));
      
      double ny_lim = sample_rate / 2.0; 
      if((abs(f1) > ny_lim) || (abs(f2) > ny_lim)) {
      }

      // We *have* our standards, you know.
      if(stopband_atten < 25.0) stopband_atten = 25.0;
  
      // set all the relevant sizes
      h.resize(num_taps);

      // create the complex image
      std::vector<std::complex<T>> H(num_taps);

      SoDa::FFT fft(h.size());
      
      // make the prototype in the frequency domain
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
      std::vector<std::complex<T>> ph(h.size());
      fft.ifft(ph, H);

      // window it
      std::vector<T> cwin(num_taps);
      SoDa::ChebyshevWindow(cwin, num_taps, stopband_atten);

      // apply gain correction
      T gain_corr = 1.0 / ((T) num_taps); 
      for(int i = 0; i < num_taps; i++) {
	ph[i] = ph[i] * cwin[i] * gain_corr;
      }

      // now do the shift
      FFT::shift(ph, ph);

      // and the reverse (makes the convolution a little faster)
      for(int i = 0; i < num_taps; i++) {
	h[i] = ph[num_taps - 1 - i];
      }
    }

    /**
     * @brief Destructor
     * 
     */
    ~FIRFilter() {
    }
    
    /** 
     * @brief Apply the filter to a single isolated buffer
     * 
     * @param out The output of the filter. 
     * @param in The input to the filter. 
     */
    void apply(std::vector<std::complex<T>> & out, std::vector<std::complex<T>> & in) {
      // we're just going to stream through the buffer one step at a time
      for(int i = 0; i < in.size(); i++) {
	Z_buffer.push(in[i]);	
	out[i] = convolveZBuf(1);
      }
    }

    // run a convolution on the Z buffer and the filter h
    // The increment is normally 1, but can be N for interpolation
    std::complex<T> convolveZBuf(int inc) {
      std::complex<T> ret(0, 0);
      for(int i = 0; i < num_taps; i += inc) {
	ret = ret + h[i] * Z_buffer[i]; 
      }
      return ret;
    }
    
    /** 
     * @brief Apply the filter to a buffer in a sequence of buffers. 
     * 
     * @param out The output of the filter. 
     * @param in The input to the filter. 
     */
    void applyCont(std::vector<std::complex<T>> & out, std::vector<std::complex<T>> & in)
    {
      apply(out, in);
    }
      

  protected:
    ///@{
    int num_taps; 

    SoDa::CircularBuffer<std::complex<T>> Z_buffer = SoDa::CircularBuffer<std::complex<T>>(0, T(0)); 

    std::vector<std::complex<T>> h;

    ///@}
    
    /**
     * @name Protected Member Functions -- shake well before using
     */
    ///@{
    
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
    ///@}
  };

}
