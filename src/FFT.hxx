#pragma once

/**
 * Simple interface to FFTW
 */
#include <complex>
#include <vector>
#include <fftw3.h>

namespace SoDa {
  class FFT {

  public:
    FFT(size_t N, unsigned int flags = FFTW_ESTIMATE | FFTW_UNALIGNED);
    
    ~FFT();

    // We could do this all with templates, but it would be even less readable that way.
    // Though the risks from code-copying are there, the simplicity of the duplicates
    // is compelling vs. the complexity of navigating all the template specializations. 
    bool fft(std::vector<std::complex<float>> & out, std::vector<std::complex<float>> & in);
    bool ifft(std::vector<std::complex<float>> & out, std::vector<std::complex<float>> & in);

    bool fft(std::vector<std::complex<double>> & out, std::vector<std::complex<double>> & in);
    bool ifft(std::vector<std::complex<double>> & out, std::vector<std::complex<double>> & in);


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
