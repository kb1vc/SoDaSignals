#pragma once
/*
  Copyright (c) 2022,2025 Matthew H. Reilly (kb1vc)
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
  Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in
  the documentation and/or other materials provided with the
  distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


///
///  @file OSFilter.hxx
///  @brief This class creates an overlap-and-save filter from a filter specification and a buffer size
///
///  This scheme replaces the awful filter creation scheme in SoDaRadio versions
///  8.x and before. 
///
///  @author M. H. Reilly (kb1vc)
///  @date   July 2022
///

#include <memory>
#include <iostream>
#include <complex>
#include <vector>
#include <fftw3.h>
#include "FilterSpec.hxx"
#include "FFT.hxx"
#include <stdexcept>
#include "Filter.hxx"

namespace SoDa {
  
  class OSFilter {
  public:
    class BadBufferSize : public std::runtime_error {
    public:
      BadBufferSize(const std::string & st, unsigned int in, unsigned int out, unsigned int req);
    };



    /**
     * constructor
     * @brief Build the filter from a filter spec for a bandpass filter
     * 
     * @param filter_spec object of class FilterSpec identifying corner frequencies and amplitudes
     * @param buffer_size the impulse response and frequency image will be this long
     * @param gain relative magnitude of input to output in the passband     
     * @param window filter window choice - we're using the window filter synthesis method. Defaults to HANN
     */
    OSFilter(FilterSpec & filter_spec, 
	     unsigned int buffer_size,
	     float gain = 1.0, 
	     Filter::WindowChoice window = Filter::HANN);        

    /**
     * @brief Alternate constructor, for very simple filters
     *
     * @param low_cutoff lower edge of the filter
     * @param high_cutoff upper edge of the filter
     * @param skirt width of transition band between cutoff and stopband
     * @param sample_rate sample rate for the input stream. 
     * @param buffer_size size of the input buffer when apply is called
     * @param stop_band_attenuation in dB
     * @param gain relative magnitude of input to output in the passband
     * @param window filter window choice - we're using the window filter synthesis method. Defaults to HANN
     */
    OSFilter(float low_cutoff, float high_cutoff, float skirt,
	     float sample_rate, unsigned int buffer_size,
	     float stop_band_attenuation = 60.0,
	     float gain = 1.0,
	     Filter::WindowChoice window = Filter::HANN);    

    

    /**
     * @brief Some subclasses of OSFilter don't have much to say
     * at construction time.
     */
    OSFilter(); 

    /**
     * @brief run the filter on a complex input stream
     * @param in_buf the input buffer I/Q samples (complex)
     * @param out_buf the output buffer I/Q samples (complex)
     * @return the length of the input buffer
     *
     * Throws OSFilter::BadRealOSFilter if the original filter spec was not "real"            
     */
    unsigned int apply(std::vector<std::complex<float>> & in_buf, 
		       std::vector<std::complex<float>> & out_buf);


    /**
     * @brief run the filter on a real input stream
     * @param in_buf the input buffer samples
     * @param out_buf the output buffer samples (this can overlap the in_buf vector)
     *
     * Throws OSFilter::BadRealOSFilter if the original filter spec was not "real"
     */
    unsigned int apply(std::vector<float> & in_buf, 
		       std::vector<float> & out_buf);


    /**
     * @brief return the number of taps in the filter.
     */
    unsigned int getTaps() { return save_buf.size() + 1; }


    /**
     * @brief return the length of the internal buffer (input buffer len + filter size - 1)
     */
    unsigned int getInternalSize() { return x_augmented.size(); }

    /**
     * @brief return what the filter object believes are its lowest and highest specified frequencies
     *
     * @return pair containing lowest, highest specified frequencies in the filter -- more or less 
     */
    std::pair<float,float> getFilterEdges() { return filter_p->getFilterEdges(); }

  private:
    void makeOSFilter(FilterSpec & filter_spec, 
		      unsigned int _buffer_size,
		      float gain, 
		      Filter::WindowChoice window);

    // This might be used in the hilbert transformer, as its filter shape
    // is rather strange. Though don't go looking for the hilbert transformer
    // here. I've never been happy with the performance around DC, and
    // its primary use was in SSB modulation, where I think the filtering
    // method probably makes much more sense and results in better
    // sideband and carrier supression. 
    void makeGenericFilter(std::vector<std::complex<float>> & H,
			   unsigned int num_taps, 
			   unsigned int buffer_size,
			   float gain = 1.0,
			   Filter::WindowChoice window_choice = Filter::HAMMING);
    
    uint32_t buffer_size; 
    
    std::unique_ptr<Filter> filter_p;
    std::vector<std::complex<float>> real_in, real_out;     
    std::vector<std::complex<float>> x_augmented;
    std::vector<std::complex<float>> save_buf;
    std::vector<std::complex<float>> y_augmented;
    std::vector<std::complex<float>> X, Y; 
  };
}

