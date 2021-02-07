#pragma once
#include <complex>
#include <vector>
#include "FFT.hxx"

namespace SoDa {
  class Filter {
  public:
    /**
     * \enum FTYPE
     * 
     * Filter type -- band pass, or band stop. 
     */
    enum FTYPE { BP, /*!< band pass */
		 BS  /*!< band stop */
    };
    
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
    Filter(FTYPE typ, int num_taps, float sample_rate, float f1, float f2, 
	   float transition_width, float stopband_atten);

    /**
     * @brief Destructor
     * 
     * Free up FFT object, if we have one. 
     */
    ~Filter() {
      if(dft_p != nullptr) {
	delete dft_p; 
      }
    }
    
    /** 
     * @brief Apply the filter to a single isolated buffer
     * 
     * @param out The output of the filter. 
     * @param in The input to the filter. 
     */
    void apply(std::vector<std::complex<float>> & out, std::vector<std::complex<float>> & in);

    /** 
     * @brief Apply the filter to a buffer in a sequence of buffers. 
     * 
     * This method implements an overlap-and-save filter, suitable for a continous 
     * signal stream. 
     * 
     * @param out The output of the filter. 
     * @param in The input to the filter. 
     */
    void applyCont(std::vector<std::complex<float>> & out, std::vector<std::complex<float>> & in);

  protected:
    ///@{
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
    ///@}
    
    /**
     * @name Protected Member Functions -- shake well before using
     */
    ///@{
    void makeImage(int image_size); 
    
    void setupNewOverlap(int invec_size);
    int save_length;
    
    void makePrototype(std::vector<std::complex<float>> & PH, float sample_rate, float f1, float f2, 
		       float transition_width); 
    
    float normalize(float freq, float sample_rate);
    ///@}
  };
}
