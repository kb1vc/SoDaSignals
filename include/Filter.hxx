#pragma once
#include <complex>
#include <vector>
#include "FFT.hxx"
#include "ChebyshevWindow.hxx"
#include <iostream>
#include <fstream>
#include <memory>
#include <list>
#include <SoDa/Format.hxx>

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
     * @param input_buffer_length If specified, the filter will be setup here. Otherwise, 
     * we'll wait until the first call to apply.
     * @param required_multiple the size of the input+overlap buffer must be a multiple of this
     */
    
    Filter(FilterType typ, int _num_taps, T sample_rate, T f1, T f2, 
	   T transition_width, T stopband_atten, 
	   int input_buffer_length = 0, 
	   int required_multiple = 1) :
      filter_type(typ), min_num_taps(_num_taps), sample_rate(sample_rate), 
      f1(f1), f2(f2), transition_width(transition_width), stopband_atten(stopband_atten), 
      last_invec_size(input_buffer_length),
      required_multiple(required_multiple)
    {
      dft_p = nullptr;
      
      if(input_buffer_length) {
	setupOverlap(last_invec_size, required_multiple);
      }
    }
    
    Filter() {
      // do nothing constructor -- used for ReSampler, f'rinstance. 
    }

    void initFilter(FilterType typ, int _num_taps, T _sample_rate, T _f1, T _f2, 
		    T _transition_width, T _stopband_atten, 
		    int _input_buffer_length = 0, 
		    int _required_multiple = 1) {

      filter_type = typ; 
      min_num_taps = _num_taps; 
      sample_rate = _sample_rate; 
      f1 = _f1; 
      f2 = _f2; 
      transition_width = _transition_width; 
      stopband_atten = _stopband_atten; 
      last_invec_size = _input_buffer_length;
      required_multiple = _required_multiple;
	
      dft_p = nullptr;
      
      if(last_invec_size) {
	setupOverlap(last_invec_size, required_multiple);
      }
      
    }

    void createFilterTaps(int num_taps) {

      // We *have* our standards, you know.
      if(stopband_atten < 25.0) stopband_atten = 25.0;

      std::cerr << "Making a prototype filter with " << num_taps << "\n";
      // set all the relevant sizes
      h.resize(num_taps);

      Proto.resize(num_taps);

      // we don't need an fft widget for the input stuff yet.
      dft_p = nullptr; 

      // make the prototype
      switch(filter_type) {
      case BP:
	makePrototype(Proto, sample_rate, f1, f2, transition_width);        
	break;
      case BS:
	makePrototype(Proto, sample_rate, f1 + transition_width, f2 - transition_width, transition_width);
	// now invert the filter
	for(auto & HE : Proto) {
	  HE = std::complex<T>(1.0, 0.0) - HE;
	}
	break;
      }

      // take the DFT
      SoDa::FFT fft(h.size());
      fft.ifft(h, Proto);
      // window it
      std::vector<T> cwin(num_taps);
      SoDa::ChebyshevWindow(cwin, num_taps, stopband_atten);

      // apply gain correction
      std::ofstream cw("chwin.out");
	
      T gain_corr = 1.0 / ((T) num_taps); 
      for(int i = 0; i < num_taps; i++) {
	cw << i << " " << h[i].real() << " " << h[i].imag() << " " << cwin[i] << " ";
	h[i] = h[i] * cwin[i] * gain_corr;
	cw << h[i].real() << " " << h[i].imag() << "\n";
      }
      cw.close();

      fft.fft(Proto, h);
      //      fft.shift(h, h);          

      dump("Proto.out", Proto);
      dump("final_h.out", h);
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
      if(H.size() != in.size()) setupOverlap(in.size(), required_multiple);
  
      // take the FFT of the input
      dft_p->fft(out, in);
      
      // now multiply by the filter image
      // and the gain correction
      T gain_corr = 1.0 / ((T) in.size());
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
      if(last_invec_size != in.size()) {
	setupOverlap(in.size(), 1); 
      }
	      
      // put the new samples at the end of the saved buffer.
      for(i = 0; i < in.size(); i++) {
	overlap_save_buffer[i + overlap_length] = in[i];
      }
      
      // take the FFT of the saved buffer
      dft_p->fft(saved_dft, overlap_save_buffer);

      dump("saved_dft_pre.out", saved_dft);      
      // multiply by the image
      for(i = 0; i < saved_dft.size(); i++) {
	auto v1 = saved_dft[i];
	saved_dft[i] = H[i] * saved_dft[i];
	writeSDFT(debug_count, i, v1, saved_dft[i], H[i]);
      }
      
      dump("saved_dft_post.out", saved_dft);

      debug_count++; 

      // take the inverse fft
      dft_p->ifft(temp, saved_dft);

      // copy the result to the output, discarding the first (overlap) results. 
      T scale = 1.0 / ((T) temp.size());
      for(i = 0; i < out.size(); i++) {
	out[i] = temp[i + overlap_length - 1] * scale; 
      }

      // save the last section of the input
      for(i = 0, j = (in.size() - overlap_length); i < overlap_length; i++, j++) {
	overlap_save_buffer[i] = in[j]; 
      }


    }

  protected:
    // set up the overlap buffer, make the filter image. 
    int findGoodSize(std::list<int> factors, int best_so_far, int test_size, int min, int max)
    {
      if(factors.empty()) return best_so_far; 
      int m = factors.front();
      std::list<int> rfactors = factors;
      rfactors.pop_front();
      
      for(int f = test_size; f < max; f *= m) {
	if((f >= min) && (f < best_so_far)) {
	  best_so_far = f;
	}
	best_so_far = findGoodSize(rfactors, best_so_far, f, min, max);
      }
      return best_so_far; 
    }
    
    void setupOverlap(int invec_size, int required_multiple)
    {
      debug_count = 0; 
  
      // how big should the transform be?
      // Make it a power of two or a power of two times 3 or 5.
      int target_min = invec_size + min_num_taps - 1;
      int target_max = target_min * 3 / 2; 

      // find a "good" size.
      // this is bone-brained stupid, but has the advantage that
      // I can explain it easily -- anything that's larger than
      // 16M entries isn't for us anyway.
      std::list<int> factors;
      factors.push_back(2); factors.push_back(3); factors.push_back(5); factors.push_back(7);
      
      int bignum = 1 << 24;
      int good_size = findGoodSize(factors, bignum, required_multiple, target_min, target_max);

      if(good_size == bignum) {
	std::stringstream ss;
	ss << "SoDa::Filter:\n\tCould not find a good overlap-save length for a " 
	   <<  min_num_taps << " tap filter and " 
	   << invec_size << "%1 entry vector\n";
	throw std::runtime_error(ss.str());
      }
  

      if((good_size - invec_size) > (invec_size / 3)) {
	std::stringstream ss;
	ss << "SoDa::Filter:\n\tGot a really bad overlap-save  length of " 
	   << good_size 
	   << " for a " << min_num_taps  << " tap filter and " 
	   << invec_size << " entry vector\n";
	throw std::runtime_error(ss.str());
      }

      saved_dft.resize(good_size);

      overlap_length = good_size - invec_size;
      overlap_save_buffer.resize(good_size);

      // clear the saved buffer
      for(int i = 0; i < overlap_save_buffer.size(); i++) {
	overlap_save_buffer[i] = std::complex<T>(0.0, 0.0);
      }

      last_invec_size = invec_size;

      int actual_num_taps = overlap_length + 1;

      createFilterTaps(actual_num_taps);
      makeImage(good_size, actual_num_taps);
    }

  protected:
	
    ///@{
    int debug_count;

    T f1, f2, transition_width, stopband_atten, sample_rate;
    int min_num_taps; 

    std::vector<std::complex<T>> h;
    std::vector<std::complex<T>> H;
    std::vector<std::complex<T>> Proto;    
    
    std::vector<std::complex<T>> overlap_save_buffer;

    std::vector<std::complex<T>> saved_dft;
    
    std::vector<std::complex<T>> temp;
    
    int last_invec_size;
    int required_multiple;

    std::shared_ptr<SoDa::FFT> dft_p;
    
    int overlap_length;

    FilterType filter_type;
    ///@}
    
    /**
     * @name Protected Member Functions -- shake well before using
     */
    ///@{
    void makeImage(int image_size, int num_taps)  {
      H.clear();

      H.resize(image_size);

      // create a new FFT object.
      dft_p = std::make_shared<SoDa::FFT>(image_size);

      std::vector<std::complex<T>> h_padded(image_size);

      int i, j;

      temp.resize(image_size);      
#if 0      
      // create the filter image.
      for(i = 0; i < image_size; i++) {
	h_padded[i] = std::complex<T>(0.0, 0.0);
      }

      for(i = 0; i < num_taps; i++) {
	h_padded[i] = h[i];
      }

      dump("h_padded.out", h_padded);
      
      // now transform it
      dft_p->fft(H, h_padded);
#else
      std::cerr << "About to make image from Proto of size " << Proto.size() << "\n";
      int lower_lim = (num_taps + 1) / 2;
      int upper_lim = image_size - (num_taps / 2);
      for(i = 0, j = 0; i < image_size; i++) {
	if(i < lower_lim) {
	  std::cerr << SoDa::Format("H[%0] = Proto[%1] = %2\n")
	    .addI(i)
	    .addI(j)
	    .addF(Proto[j].real());
	  H[i] = Proto[j];
	  j++; 
	}
	else if(i > upper_lim) {
	  std::cerr << SoDa::Format("H[%0] = Proto[%1] = %2\n")
	    .addI(i)
	    .addI(j)
	    .addF(Proto[j].real());
	  H[i] = Proto[j];
	  j++; 
	}
	else {
	  H[i] = std::complex<T>(0.0, 0.0);
	}
      }
      std::cerr << "Made image\n";      
#endif      
      dump("H.out", H);
    }
    
    void makePrototype(std::vector<std::complex<T>> & PH, T sample_rate, T f1, T f2, 
		       T transition_width) {
      // Assume we're building a band pass filter
      T lo_stop = f1 - transition_width;
      T lo_pass = f1;
      T hi_pass = f2;
      T hi_stop = f2 + transition_width;

      int num_taps = PH.size();
      T bucket_width = sample_rate / ((T) num_taps);

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

      std::cerr << SoDa::Format("f1 %0 f2 %1 lo_stop %2 lo_pass %3 hi_pass %4 hi_stop %5\n")
	.addF(f1)
	.addF(f2)
	.addF(lo_stop)
	.addF(lo_pass)
	.addF(hi_pass)
	.addF(hi_stop);
      
      for(int i = 0; i < num_taps; i++) {
	int f_idx = i - (num_taps / 2);
	T freq = bucket_width * ((T) f_idx); 
	T v; 
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

      dump("PH.out", PH);
      FFT::shift(PH, PH);
      dump("PH_S.out", PH);
    }

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

    
    T normalize(T freq, T sample_rate) {
      T res = freq / sample_rate;

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
