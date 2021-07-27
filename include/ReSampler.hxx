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
      out_buffer_len = (in_buffer_len * interpolate) / decimate;

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
      lpf.initFilter(SoDa::FilterType::BP, 41, 
		     sample_rate, -freq_corner, freq_corner,
		     0.1 / interpolate, // transition_width * sample_rate,
		     50, 
		     in_buffer_len, 
		     decimate);
      
      
      
      int interpolate_buffer_size = in_buffer_len * interpolate;

      out_tmp.resize(interpolate_buffer_size);
      int_buffer.resize(interpolate_buffer_size);

      // zero out the interpolated vector
      for(int i = 0; i < interpolate_buffer_size; i++) {
	int_buffer[i] = std::complex<T>(0.0, 0.0);
      }
      
      debug_count = 0; 
    }

    /**
     * @brief Destructor
     * 
     */
    ~ReSampler() {
    }
    
    /** 
     * @brief Apply the Resampler to a single isolated buffer : we don't do that.
     * 
     * @param out The output of the filter. 
     * @param in The input to the filter. 
     */
    void apply(std::vector<std::complex<T>> & out, std::vector<std::complex<T>> & in) {
      std::ofstream tmp("rsinc_in.dat");
      for(int i = 0; i < in.size(); i++) {
	tmp << in[i].real() << " " << in[i].imag() << "\n";
      }
      tmp.close();
      applyCont(out, in);
    }

    void findMax(std::vector<std::complex<T>> & v, const std::string & name) {
      int imax = 0;
      T vmax = 0.0;
      for(int i = 0; i < v.size(); i++) {
	T ma = abs(v[i]);
	if(ma > vmax) {
	  vmax = ma;
	  imax = i; 
	}
      }
      std::cerr << SoDa::Format("Max at %0[%1] = %2\n")
	.addS(name)
	.addI(imax)
	.addF(vmax);
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
      std::cerr << "In applyCont\n";
      
      if(in.size() != in_buffer_len) {
	throw SoDa::ReSampler_MismatchBufferSize("In", in_buffer_len, in.size());
      }
      if(out.size() != out_buffer_len) {
	throw SoDa::ReSampler_MismatchBufferSize("Out", out_buffer_len, out.size());	
      }

      // interpolate
      for(i = 0, j = 0; i < in.size(); i++, j += interpolate) {
	int_buffer[j] = in[i];
      }

      dump("invec.out", in);
      dump("int_buffer.out", int_buffer);
      
      std::cerr << "Applying filter\n";
      // apply the filter.
      lpf.applyCont(out_tmp, int_buffer);
      std::cerr << "Applied filter\n";

      dump("out_temp.out", out_tmp);
      
      // now stride through and downsample ?
      for(i = 0; i < out.size(); i++) {
	out[i] = out_tmp[i * decimate];
      }

      //      lpf.applyCont(out, out);
      
      dump("out.out", out);
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
    
    std::vector<std::complex<T>> int_buffer;
    
    int out_buffer_len, in_buffer_len; 

    int interpolate, decimate; 
    
    std::vector<std::complex<T>> out_tmp;

    SoDa::Filter<T> lpf; // lowpass filter for interp/decimation
    ///@}
  };

}
