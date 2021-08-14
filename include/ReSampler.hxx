#pragma once
#include <complex>
#include <vector>
#include "Filter.hxx"
#include <iostream>
#include <memory>
#include <stdexcept>
#include <SoDa/Format.hxx>
#include <iostream>
#include <fstream>

namespace SoDa {
  
  class ReSampler_BadBufferSize : public std::runtime_error {
  public:
    ReSampler_BadBufferSize(int in_blen, int interp, int dec) :
      std::runtime_error(SoDa::Format("SoDa::ReSampler Buffer Length %0 can't be rationally resampled as %1/%2\n")
			 .addI(in_blen)
			 .addI(interp)
			 .addI(dec)
			 .str()) {
    }
  };

  class ReSampler_MismatchBufferSize : std::runtime_error {
  public:
    ReSampler_MismatchBufferSize(const std::string dir, 
				 int in_blen, 
				 int exp_in_blen) :
      std::runtime_error(SoDa::Format("SoDa::ReSampler %0put Buffer Length %1 should have been %2\n")
			 .addS(dir)
			 .addI(in_blen)
			 .addI(exp_in_blen)
			 .str()) {
    }
  };
  
  template <typename T>
  class ReSampler {
  public:
    
    /**
     * @brief Constructor -- A Rational Resampler in the frequency domain
     * 
     * This builds on the overlap-and-save SoDa::Filter
     * 
     * @param in_buffer_len length of input buffer to be resampled
     * @param interpolate interpolation (upsampling) factor
     * @param decimate decimation (downsampling) factor
     * @param transition_width width of skirt between output low-pass freq
     * and output nyquist rate
     *
     * This constructor will throw a "SoDa::Resampler_BadRate" exception
     * if the up/down ratio can't be implemented in the specified buffer
     * lengths. 
     *
     * ## The Algorithm
     *
     * This is a fairly straightforward filter-based implementation of a 
     * rational resampler. It takes place in two stages: Upsampling, 
     * followed by Dowsampling. Most of the heavy lifting is in the frequency
     * domain.  Interpolate by I, decimate by D
     * 
     * -# x'[nI] = x[n] -- other X'[m] = 0
     * -# y' = filter(x')
     * -# y[n] = y[n D]
     *
     * ## The implementation
     * 
     * We'll use an overlap-and-save method, building on the components of the
     * SoDa::Filter class. The LPF filter H[] will be created in the Filter. Its
     * bandwidth will be set by the Interpolate and Decimate factors so that the
     * aliasing in both steps is addressed.  The cutoff frequency needs to be
     * low enough to provide the interpolation smoothing for the upsampling and
     * the alias elimination for the downsampling. 
     * $F_{lpf} = F_{s} I} / (2 D)$ symetric about 0. (really a BPF)
     */
    ReSampler(int in_buffer_len, 
	      int interpolate,
	      int decimate,
	      float transition_width) : 
      in_buffer_len(in_buffer_len), interpolate(interpolate), decimate(decimate)
    {
      transition_width = 0.1;
      int out_buffer_len = (in_buffer_len * interpolate) / decimate;

      if(((out_buffer_len * decimate) / interpolate) != in_buffer_len) {
	// then we screwed up and the ratio won't work. 
	// yes, this is overly restrictive, but why buy agony wholesale?
	throw SoDa::ReSampler_BadBufferSize(in_buffer_len, interpolate, decimate); 
      }

      // we need a filter.  Nothing fancy.  An LPF at 0.5 / I - TW
      // sample rate is arbitrary, make it 1M
      float sample_rate = 1;
      float freq_corner = 0.4 / ((float) interpolate); // (0.4 /((float) interpolate) - transition_width);
      

      // now setup the image and overlap-and-save buffers
      int overlap_mul = decimate;
      // we need an even overlap to make the zero stuffing
      // symmetric
      if(overlap_mul & 1) { // if it is odd
	overlap_mul = overlap_mul * 2; 
      }

      int diff_len = out_buffer_len - in_buffer_len;
      
      lpf.initFilter(SoDa::FilterType::BP, 41, 
		     sample_rate, -freq_corner, freq_corner,
		     0.1 * freq_corner,
		     50, 
		     in_buffer_len,
		     // overlap buffer length must be an integer multiple of
		     // interp / decimate
		     [diff_len, overlap_mul] (int s) ->
		     bool { return ( ( s - diff_len) % overlap_mul) == 0; });

      
      
      
      int pre_image_size = in_buffer_len + lpf.getOverlapLen();
      pre_image.resize(pre_image_size); 

      int interp_image_size = out_buffer_len + lpf.getOverlapLen() * interpolate / decimate;
      interp_image.resize(interp_image_size); 

      resamp_res.resize(interp_image_size);

      dft_p = new SoDa::FFT(interp_image_size);
      
      debug_count = 0; 
    }

    /**
     * @brief Destructor
     * 
     */
    ~ReSampler() {
    }
    
    /** 
     * @brief Apply the Resampler to a buffer stream
     * 
     * @param out The output of the filter. 
     * @param in The input to the filter. 
     */
    void applyCont(std::vector<std::complex<T>> & out, std::vector<std::complex<T>> & in) {
      std::ofstream tmp("rsinc_in.dat");

      // do the anti-aliasing filter
      lpf.applyContFFT(pre_image, in);

      // how many zeros do we need to stuff? 
      auto M = lpf.getOverlapLen();
      int discard = M * interpolate / decimate;      
      int stufflen = discard + out.size() - in.size();
      int half = (in.size() + M) / 2;
      // now the zero stuffing of M * interp / decimation + (out.len - in.len) zeros
      int j = 0;
      for(int i = 0; i < interp_image.size(); i++) {
	if(i < half) {
	  interp_image[i] = pre_image[j++];
	}
	else if(i < (half + stufflen)) {
	  interp_image[i] = std::complex<T>(0.0, 0.0);
	}
	else {
	  interp_image[i] = pre_image[j++];	  
	}
      }

      // now ifft
      dft_p->ifft(resamp_res, interp_image);


      // discard first M * interp / decimation samples.
      for(int i = 0; i < out.size(); i++) {
	out[i] = resamp_res[i + discard];
      }

      // return. 
    }

      

  protected:

    void dump(const std::string & fn, std::vector<std::complex<T>> vec) {
      std::ofstream out(fn);

      for(auto & v : vec) {
	T mag = v.real() * v.real() + v.imag() * v.imag();
	out << v.real() << " " << v.imag() << " " << mag << "\n";
      }

      out.close();
    }
    
    ///@{
    int debug_count;
    
    int interpolate, decimate; 
    int in_buffer_len; 

    std::vector<std::complex<T>> pre_image;
    std::vector<std::complex<T>> interp_image;
    std::vector<std::complex<T>> resamp_res;
    
    
    SoDa::Filter<T> lpf; // lowpass filter for interp/decimation
    SoDa::FFT * dft_p;
    
    
    ///@}
  };

}
