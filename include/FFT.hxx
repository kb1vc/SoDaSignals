#pragma once
#include <complex>
#include <vector>
#include <fftw3.h>
/**
 * @file FFT.hxx
 * 
 * @author M. H. Reilly (kb1vc)
 * @date Jan2021
 */
namespace SoDa {
/**
 * @class SoDa::FFT
 * 
 * @brief Simple interface to FFTW -- provides DFT and inverse DFT functions
 *
 */
  class FFT {
  public:
    /**
     * @brief constructor
     * 
     * Build an object that can produce DFT/IDFT for complex vectors 
     * of a particular size.  By constructing an FFT object to be applied
     * to multiple vectors, we save the initialization and measurement
     * costs associated with FFTW's run-time optimization. 
     * 
     * @param N the number of elements in an input vector passed to fft or ifft
     * @param flags FFTW optimization flags.  FFTW_UNALIGNED is, unless the user
     * knows better, all-but-mandatory.  FFTW_ESTIMATE provides reasonable performance
     * for modern X86 processors, and is sufficiently close to the result from 
     * FFTW_MEASURE as to be "just fine."
     * 
     */
    FFT(size_t N, unsigned int flags = FFTW_ESTIMATE | FFTW_UNALIGNED);
    
    ~FFT();

    // We could do this all with templates, but it would be even less readable that way.
    // Though the risks from code-copying are there, the simplicity of the duplicates
    // is compelling vs. the complexity of navigating all the template specializations. 
  
    /**
     * @brief forward DFT on vector of complex<float>
     * 
     * @param out the output vector.  It must have exactly as many elements as 
     * were passed as the N parameter of the constructor.
     * @param in the input vector.  It must have exactly as many elements as 
     * were passed as the N parameter of the constructor.
     * 
     * @return 'true' all the time. 
     * 
     * @throw runtime_exception If the input or output vector is not equal to the
     * "N" parameter passed in to the constructor, any fft or ifft function will
     * throw an exception. 
     */
    bool fft(std::vector<std::complex<float>> & out, std::vector<std::complex<float>> & in);

    /**
     * @brief inverse DFT on vector of complex<float>
     * 
     * @param out the output vector.  It must have exactly as many elements as 
     * were passed as the N parameter of the constructor. 
     * @param in the input vector.  It must have exactly as many elements as 
     * were passed as the N parameter of the constructor.
     * 
     * @return 'true' all the time. 
     * 
     * @throw runtime_exception If the input or output vector is not equal to the
     * "N" parameter passed in to the constructor, any fft or ifft function will
     * throw an exception. 
     */
    bool ifft(std::vector<std::complex<float>> & out, std::vector<std::complex<float>> & in);

    /**
     * @brief forward DFT on vector of complex<double>
     * 
     * @param out the output vector.  It must have exactly as many elements as 
     * were passed as the N parameter of the constructor.
     * @param in the input vector.  It must have exactly as many elements as 
     * were passed as the N parameter of the constructor.
     * 
     * @return 'true' all the time. 
     * 
     * @throw runtime_exception If the input or output vector is not equal to the
     * "N" parameter passed in to the constructor, any fft or ifft function will
     * throw an exception. 
     */
    bool fft(std::vector<std::complex<double>> & out, std::vector<std::complex<double>> & in);
    /**
     * @brief inverse DFT on vector of complex<double>
     * 
     * @param out the output vector.  It must have exactly as many elements as 
     * were passed as the N parameter of the constructor.
     * @param in the input vector.  It must have exactly as many elements as 
     * were passed as the N parameter of the constructor.
     * 
     * @return 'true' all the time. 
     * 
     * @throw runtime_exception If the input or output vector is not equal to the
     * "N" parameter passed in to the constructor, any fft or ifft function will
     * throw an exception. 
     */
    bool ifft(std::vector<std::complex<double>> & out, std::vector<std::complex<double>> & in);


    /**
     * @brief static utility function to "rotate" an FFT result 
     * 
     * DFT operations will return a vector with the DC component of the transform
     * at index 0, the nyquist frequency at N/2 and the lowest negative component
     * at index N-1. 
     *
     * Earthlings like looking at spectrum displays with the DC component in the middle, 
     * the negative frequencies to the left, and the positive frequencies on the
     * right. This static function does that.  Note that after this transformation, 
     * the result of an ifft or fft operation on the resulting vector may
     * be a surprise. 
     *
     * @param out the output vector.
     * @param in the input vector. 
     * 
     */
    template<typename T> static void shift(std::vector<T> & out, const std::vector<T> & in) {
      // make this work even if out and in are the same. 
      unsigned int N = in.size();
      if(N % 2 == 0) {
	// easy, swap
	int i, j; 
	for(i = 0, j = N/2; i < N / 2; i++, j++) {
	  T v = in[j];
	  out[j] = in[i];
	  out[i] = v; 
	}
      }
      else {
	// this is an odd one... ;)
	T last_v = in[0];
	int bump = (N - 1) / 2;
	int cur_idx = bump;
	for(int i = 0; i < N; i++) {
	  T v = in[cur_idx];
	  out[cur_idx] = last_v;
	  last_v = v;
	  cur_idx = (cur_idx + bump) % N;
	}
      }
    }

  private:
    /**
     * @brief throw an exception if either the input or output buffer size does
     * not match the plan
     */
    void checkInOut(size_t outsize, size_t insize);
    
    fftwf_plan fplan_f;
    fftwf_plan fplan_i;    

    fftw_plan dplan_f;
    fftw_plan dplan_i;    
    
    // use these to create the basic plans.  sigh. 
    fftwf_complex * f_dummy_in, * f_dummy_out;    
    fftw_complex * d_dummy_in, * d_dummy_out;
    
    size_t dim;
  };
}