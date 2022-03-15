#pragma once
#include <complex>
#include <vector>
#include <functional>
#include "FFT.hxx"
#include "Windows.hxx"
#include <iostream>
#include <fstream>
#include <memory>
#include <list>
#include <SoDa/Format.hxx>

/**
 * Generic Filter Class
 * 
 * Build a FIR overlap-and-save filter from a passband/stopband specification over the range 
 * - fsamp / 2 to fsamp / 2
 * 
 * The process: 
 * 
 * 1. Given a minimum filter tap size of N, and an input buffer size of B find a 
 * convenient FFT length Q that is greater than B + N.  The overlap size is 
 * M = Q - N
 * 2. Build a frequency domain image of the filter Hp for an M + 1 tap filter over the
 * range -Fs/2 to Fs/2. 
 * 2a. ifft to make hp
 * 2b. window hp with an M + 1 element hamming window. 
 * 2c. fft to make Hwp
 * 3. Expand Hwp to H -- a Q element filter by stuffing zeros. 
 * 
 * The applyCont method performs the overlap-and-save calculation using the filter H. 
 */
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
     * @param sample_rate Samples per second  Fs
     * @param f1 The lower pass/stop band edge in Hz    -0.5 Fs < f1 < 0.5 Fs
     * @param f2 The upper pass/stop band edge in Hz    -0.5 Fs < f2 < 0.5 Fs  f1 < f2
     * @param transition_width  Width of skirts (more or less)
     * @param input_buffer_length we fix the input chunk length here. 
     */
    
    Filter(FilterType typ, int _num_taps, T sample_rate, T f1, T f2, 
	   T transition_width, T stopband_atten,
	   int input_buffer_length) : low_freq(f1), high_freq(f2)
    {
      initFilter(typ, _num_taps, 
		 sample_rate, 
		 f1 / sample_rate, 
		 f2 / sample_rate, 
		 transition_width / sample_rate, 
		 input_buffer_length); 
    }
    
    Filter() {
      // do nothing constructor -- used for ReSampler, f'rinstance.
    }

    void initFilter(FilterType typ, int min_num_taps, T sample_rate, T f_lo, T f_hi, 
		    T transition_width,
		    int buffer_len, 
		    std::function<bool(const int)> test = 
		    [] (int s) -> bool { return true; })
    {

      debug_count = 0; 
      // Find a good size:  This is an FFT size > buffer_len + num_taps
      int target = buffer_len + min_num_taps; 
      int pow_of_two = 1;
      while(pow_of_two < target) {
	pow_of_two *= 2; 
      }

      // now do the factor search
      int fft_len = findGoodSize(std::list<int> { 2, 3, 5, 7 }, pow_of_two, target, 1, test);

      // Now we know the size of the overlap.
      overlap_len = fft_len - buffer_len; 

      // create a frequency domain image of the filter
      int num_taps = overlap_len + 1;
      auto Hi = createImage(num_taps, sample_rate, f_lo, f_hi, transition_width, typ);
      
      // now apply a chebyshev window to it
      auto hwp = windowImpulse(Hi);
      
      // expand the filter
      H = expandFilter(hwp, fft_len);

      // create the "save buffer" -- this is the buffer that
      // will hold the overlap + the new buffer. 
      overlap_save_buffer.resize(fft_len); 

      // create the FFT object
      dft_p = new SoDa::FFT(fft_len); 
      
      // setup buffers for the intermediate values
      saved_dft.resize(fft_len); 
      saved_idft.resize(fft_len); 
      
      // make sure the amplitude is preserved
      scale_factor = 1.0 / T(fft_len);

    }

    
    /**
     * @brief create a frequency domain image of a filter
     * 
     * All frequencies are normalized to -fs to fs
     */
    
    std::vector<std::complex<T>> createImage(int num_taps, 
					     T samp_rate, 
					     T f_lo,
					     T f_hi,
					     T trans,
					     FilterType typ) 
    {
      std::vector<std::complex<T>> H(num_taps);

      // assume a band-pass filter
      T lo_stop = f_lo - 0.5 * trans;
      T lo_pass = f_lo;
      T hi_pass = f_hi;
      T hi_stop = f_hi + 0.5 * trans;

      T f_incr = 1.0 / num_taps;

      int i;
      double f; 
      for(i = 0, f = -0.5 ; i < num_taps; i++, f += f_incr) {
	if((f < lo_stop) || (f > hi_stop)) {
	  H[i] = std::complex<T>(0.0, 0.0);
	}
	else if((f > lo_pass) && (f < hi_pass)) {
	  H[i] = std::complex<T>(1.0, 0.0);
	}
	else if(f < lo_pass) {
	  // transition band.
	  H[i] = std::complex<T>((f - lo_stop) * trans, 0.0);
	}
	else if(f > hi_pass) {
	  // transition band.
	  H[i] = std::complex<T>(1.0 - (f - hi_stop) * trans, 0.0);
	}
      }


      if(typ == FilterType::BS) {
	// this is a band-stop filter.  Invert the image
	auto one = std::complex<T>(1.0, 0.0);
	for(int i = 0; i < num_taps; i++) {
	  H[i] = one - H[i]; 
	}
      }

      // now do an FFT shift
      SoDa::FFT fft(num_taps);
      fft.shift(H, H);
      
      return H; 
    }

    std::vector<std::complex<T>> windowImpulse(std::vector<std::complex<T>> & Hi) {
      int flen = Hi.size();

      std::vector<std::complex<T>> hi(flen);

      // first make the impulse response
      SoDa::FFT fft(flen);
      fft.ifft(hi, Hi);

      // now create a hamming window
      std::vector<T> cwin(flen);
      SoDa::HammingWindow(cwin, flen);

      // apply the window
      T sf = 1.0 / T(flen);
      for(int i = 0; i < flen; i++) {
	hi[i] = hi[i] * cwin[i];
      }

      fft.shift(hi, hi);
      
      return hi;
    }

    /**
     * Take a prototype windowed filter Hp and expand it to a filter of length [new_len] by zero padding. 
     *
     * @param hi the impulse response of the prototype (short) filter
     * @param new_len the length of the new filter (in an OAS filter this is "Q" long)
     * @return H a new_len tap filter (in the frequency domain) 
     */
    std::vector<std::complex<T>> expandFilter(std::vector<std::complex<T>> & hi, int new_len) {
      int short_len = hi.size();

      // now copy the impulse response to a bigger h
      std::vector<std::complex<T>> h(new_len);
      // zero out the impulse response
      for(int i = 0; i < new_len; i++) {
	h[i] = std::complex<T>(0.0, 0.0); 	
      }
      // now fill
      for(int i = 0; i < short_len; i++) {
	int j = i;
	h[j] = hi[i]; 
      }
      
      // now expand.
      std::vector<std::complex<T>> H(new_len);
      SoDa::FFT fft(new_len);
      fft.fft(H, h);

      // scale it
      T scale = 1.0 / (T(short_len) * T(new_len));
      for(auto & Hk : H) {
	Hk = Hk * scale; 
      }
      
      return H; 
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
     * @brief Apply the filter to a buffer in a sequence of buffers, but don't 
     * apply IFFT to the frequency domain result. 
     * 
     * This method implements an overlap-and-save filter, suitable for a continous 
     * signal stream. 
     * 
     * @param out The frequency domain output of the filter. 
     * @param in The input to the filter. 
     */
    void applyContFFT(std::vector<std::complex<T>> & out, std::vector<std::complex<T>> & in) {
      int i, j;

      // put the new samples at the end of the saved buffer.
      for(i = 0; i < in.size(); i++) {
	overlap_save_buffer.at(i + overlap_len) = in[i];
      }

      // take the FFT of the saved buffer
      dft_p->fft(out, overlap_save_buffer);

      // multiply by the image
      for(i = 0; i < out.size(); i++) {
	out[i] =  H[i] * out[i];
      }
      
      debug_count++; 


      // save the last section of the input
      for(i = 0, j = (in.size() - overlap_len); i < overlap_len; i++, j++) {
	overlap_save_buffer[i] = in[j]; 
      }
    }

    /** 
     * @brief Apply the filter to a buffer in a sequence of buffers. 
     * 
     * This method implements an overlap-and-save filter, suitable for a continous 
     * signal stream. 
     * 
     * @param out The output of the filter. 
     * @param in The input to the filter. 
     * @param decimation normally 1, but otherwise skip this many elements 
     * in the idft buffer between saved samples in out.  Used for resamplers. 
     */
    void applyCont(std::vector<std::complex<T>> & out, 
		   std::vector<std::complex<T>> & in, 
		   int decimation = 1)
    {
      applyContFFT(saved_dft, in); 
      
      // take the inverse fft
      dft_p->ifft(saved_idft, saved_dft);

      // copy the result to the output, discarding the first (overlap) results. 
      for(int i = 0; i < out.size(); i++) {
	out[i] = saved_idft.at(i * decimation + overlap_len);
      }
    }

    std::pair<T,T> getFilterEdges() {
      return std::pair<T,T>(low_freq, high_freq); 
    }
    
    unsigned int getOverlapLen() { return overlap_len; }
    
    unsigned int getSaveBufLen() { return saved_idft.size(); }
  protected:
    /**
     * @brief find a good FFT size. 
     *
     * Step 1... given a minimum filter tap size of N, and an input buffer of B find a
     *  convenient FFT length Q that is greater than B + N, but not too much so.
     *
     * This is a recursive function: given a list of factors that are useful in mixed-radix
     * FFTs, look for something good.  The result should be the smallest Q that is a product
     * of one or more of the supplied factors
     * 
     * @param factors a list of factors for Q 
     * @param best_so_far the best candidate Q we've seen in the search so far.  start out with 2(B+N) we'll do better.
     * @param target (N + B)
     * @param C current length -- not useful until C > target
     * @param isAcceptable function returning true iff the calculated Q meets some criteria (e.g. it is a multiple of 9)
     * @return best length so far.  we don't need the factors, just the length
     */
    int findGoodSize(std::list<int> factors, int best_so_far, int target, int C,
		     std::function<bool(int)> isAcceptable)
    {

      if(factors.empty()) return best_so_far; 
      int m = factors.front();
      std::list<int> rfactors = factors; // remaining factors
      rfactors.pop_front();

      for(int f = C; f < best_so_far; f *= m) {
	if((f >= target) && (f < best_so_far) && isAcceptable(f)) {
	  best_so_far = f;
	}
	best_so_far = findGoodSize(rfactors, best_so_far, target, f, isAcceptable);
      }
      return best_so_far; 
    }
    

  protected:
	
    int debug_count;

    std::vector<std::complex<T>> H;

    int overlap_len; 

    std::vector<std::complex<T>> overlap_save_buffer;
    std::vector<std::complex<T>> saved_dft;
    std::vector<std::complex<T>> saved_idft;    

    T scale_factor; 

    SoDa::FFT * dft_p;

    T low_freq, high_freq;
    
    void writeSDFT(int count, int idx, std::complex<T> & vin, std::complex<T> & vout, std::complex<T> & HI) {
      auto min = vin.real() * vin.real() + vin.imag() * vin.imag();
      auto mout = vout.real() * vout.real() + vout.imag() * vout.imag();
      auto hmag = HI.real() * HI.real() + HI.imag() * HI.imag();
      std::cout << "JJ " << count << " " << idx << " " << min << " " << mout << " " << hmag << "\n";
    }


    void dump(const std::string & fn, std::vector<std::complex<T>> vec) {
      std::ofstream out(fn);

      for(auto & v : vec) {
	T mag = v.real() * v.real() + v.imag() * v.imag();
	out << v.real() << " " << v.imag() << " " << mag << "\n";
      }

      out.close();
    }


    void dump(const std::string & fn, std::vector<T> vec) {
      std::ofstream out(fn);

      for(auto & v : vec) {
	out << v << "\n";
      }

      out.close();
    }
    
    
    T normalize(T freq, T sample_rate) {
      T res = freq / sample_rate;

      res = M_PI * res;

      if(res > M_PI) res = res - 2.0 * M_PI;
      if(res < -M_PI) res = res + 2.0 * M_PI; 
      return res; 
    }
  };

  typedef Filter<float> FilterF;
  typedef Filter<double> FilterD;
}
