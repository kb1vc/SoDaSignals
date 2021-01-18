#pragma once
#include <complex>
#include <vector>
#include "FFT.hxx"

namespace SoDa {
  class Filter {
  public:
    enum FTYPE { LO, HI, BAND };
    
    Filter(FTYPE typ, int num_taps, float sample_rate, float f1, float f2, 
	   float transition_width, float stopband_atten);

    void getTaps(std::vector<float> & h);
    
    void getTransform(std::vector<float> & H);

    void apply(std::vector<std::complex<float>> & out, std::vector<std::complex<float>> & in);

    void applyCont(std::vector<std::complex<float>> & out, std::vector<std::complex<float>> & in);

    template<typename T> 
    friend void priv_apply(std::vector<T> & out, std::vector<T> & in, 
						SoDa::Filter & f);
    template<typename T, typename T2> 
    friend void priv_applyCont(std::vector<T> & out, std::vector<T> & in,
			       std::vector<T2> & saved, std::vector<T2> saved_dft,
			       SoDa::Filter & f);
  protected:
    int debug_count;
    
    int num_taps; 

    std::vector<std::complex<float>> h;
    std::vector<std::complex<float>> H;
    
    std::vector<std::complex<float>> saved;

    std::vector<std::complex<float>> saved_dft;
    
    std::vector<std::complex<float>> temp;
    
    int last_invec_size;

    float bucket_width;

    SoDa::FFT * dft_p; 
    
    void makeImage(int image_size); 
    
    void setupNewOverlap(int invec_size);
    int save_length;
    
    void makePrototype(std::vector<std::complex<float>> & PH, float sample_rate, float f1, float f2, 
		       float transition_width); 
    
    float normalize(float freq, float sample_rate);
  };
}
