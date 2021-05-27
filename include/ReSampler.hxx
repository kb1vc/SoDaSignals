#pragma once
#include <complex>
#include <vector>
#include "FFT.hxx"
#include "Filter.hxx"
#include <iostream>
#include <memory>
#include <stdexcept>
#include <SoDa/Format.hxx>

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
  class ReSampler :  Filter<T> {
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
     * The key here is to pick an overlap-and-save size that ends up a multiple of the 
     * output buffer length.  Otherwise resampling artifacts can kick in and cause 
     * discontinuities from one buffer to the next. 
     */
    ReSampler(int in_buffer_len, 
	      int interpolate,
	      int decimate,
	      float transition_width) : 
      Filter<T>(), in_buffer_len(in_buffer_len), interpolate(interpolate), decimate(decimate)
    {
      out_buffer_len = (in_buffer_len * interpolate) / decimate;

      if(((out_buffer_len * decimate) / interpolate) != in_buffer_len) {
	// then we screwed up and the ratio won't work. 
	// yes, this is overly restrictive, but why buy agony wholesale?
	throw SoDa::ReSampler_BadBufferSize(in_buffer_len, interpolate, decimate); 
      }

      // we need a filter.  Nothing fancy.  An LPF at 0.5 * I/D - TW
      // sample rate is arbitrary, make it 1M
      float sample_rate = 1e6; 
      float freq_corner = sample_rate * 
	(0.5 /((float) decimate) - transition_width);

      // Keep the number of taps to about 128 -- use Fred Harris' rule of thumb
      int num_taps = static_cast<int>(round(60.0 / (22.0 * transition_width)));
      
      // we need the filter stuff to create the LPF/BPF image
      // and the overlap-and-save buffers. 
      Filter<T>::initFilter(FilterType::BP, 
			    num_taps, 
			    sample_rate, 
			    -freq_corner, freq_corner, 
			    transition_width, 50);

      // now setup the image and overlap-and-save buffers
      Filter<T>::setupOverlap(in_buffer_len * interpolate, decimate); 

      // create the FFT widgets
      std::cerr << SoDa::Format("out_buffer_len %0  saved_buffer_len %1 in_buffer_len %2\n")
	.addI(out_buffer_len)
	.addI(Filter<T>::saved_dft.size())
	.addI(in_buffer_len);
      
      int saved_buffer_size = Filter<T>::saved_dft.size();
      dft_interp_p = std::make_shared<SoDa::FFT>(saved_buffer_size);
      dft_dec_p = std::make_shared<SoDa::FFT>(saved_buffer_size);

      out_dft.resize(saved_buffer_size);
      out_tmp.resize(saved_buffer_size);
      
      // zero out the saved vector
      for(int i = 0; i < Filter<T>::saved.size(); i++) {
	Filter<T>::saved[i] = std::complex<T>(0.0, 0.0);
      }
      
      debug_count = 0; 
    }

    /**
     * @brief Destructor
     * 
     * Free up FFT object, if we have one. 
     */
    ~ReSampler() {
    }
    
    /** 
     * @brief Apply the Resampler to a single isolated buffer : we don't do that..
     * 
     * @param out The output of the filter. 
     * @param in The input to the filter. 
     */
    void apply(std::vector<std::complex<T>> & out, std::vector<std::complex<T>> & in) {
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

      if(out.size() != out_buffer_len) {
	throw SoDa::ReSampler_MismatchBufferSize("In", in_buffer_len, in.size());
      }
      if(out.size() != out_buffer_len) {
	throw SoDa::ReSampler_MismatchBufferSize("Out", out_buffer_len, out.size());	
      }

      
      // put the new samples at the end of the saved buffer.
      // stuffing them as we go along. 
      for(i = 0; i < in.size(); i++) {
	Filter<T>::saved[Filter<T>::save_length + i * interpolate] = in[i];
      }

      std::cerr << "saved_dft.size() = " << Filter<T>::saved_dft.size() << "\n";
      
      // take the FFT of the saved buffer
      std::cerr << "interpolating\n";
      dft_interp_p->fft(Filter<T>::saved_dft, Filter<T>::saved);

      findMax(Filter<T>::saved_dft, "saved_dft");
      
      for(int i = Filter<T>::save_length; i < Filter<T>::saved.size(); i++) {
	std::cout << "ZZ " << debug_count++ << " " << Filter<T>::saved[i].real() << "\n";
      }

      // multiply by the image
      // but skip the parts we don't care about
      int lim_low = out_buffer_len / 2;
      int lim_high =  Filter<T>::saved_dft.size() - (out_buffer_len / 2); // is this right?
      int si = 0;

      for(i = 0; i < Filter<T>::saved_dft.size(); i++) {
	out_dft[i] = Filter<T>::H[i] * Filter<T>::saved_dft[i];	
      }

      findMax(out_dft, "out_dft");
      
      // take the inverse fft
      std::cerr << "decimating\n";
      dft_dec_p->ifft(out_tmp, out_dft);
      // now stride through and downsample ?
      for(i = 0; i < out.size(); i++) {
	out[i] = out_tmp[Filter<T>::save_length + i * decimate];
      }

      std::cerr << "Saving  save length = " << Filter<T>::save_length << "\n";
      // save the last section of the input
      // stuffing it as we go along
      for(i = 0, j = (in.size()- (Filter<T>::save_length / interpolate)); 
	  (i * interpolate) < Filter<T>::save_length; i++, j++) {
	Filter<T>::saved[i * interpolate] = in[j]; 
      }
    }
      

  protected:
    ///@{
    int debug_count;
    
    std::vector<std::complex<T>> temp;
    
    int out_buffer_len, in_buffer_len; 

    int interpolate, decimate; 
    
    std::shared_ptr<SoDa::FFT> dft_interp_p, dft_dec_p;

    std::vector<std::complex<T>> out_dft, out_tmp; 
    ///@}
  };

}
