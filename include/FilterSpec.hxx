#pragma once
/*
  Copyright (c) 2022, 2025 Matthew H. Reilly (kb1vc)
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
///  @file FilterSpec.hxx
///  @brief This class describes a filter profile -- corner frequencies and
///  amplitudes. 
///
///  This scheme replaces the awful filter creation scheme in SoDaRadio versions
///  8.x and before. 
///
///  @author M. H. Reilly (kb1vc)
///  @date   July 2022
///

#include <vector>
#include <list>
#include <complex>
#include <stdexcept>
#include <SoDa/Format.hxx>

namespace SoDa {

  /**
   * @class FilterSpec
   *
   * @brief Describes the shape of a filter. This is used by SoDa::Filter and SoDa::OSFilter
   * to describe complicated filters that aren't band-pass, band-stop, low-pass, or high-pass. 
   */
  class FilterSpec {
  public:
    
    enum FType { 
      REAL, ///< The filter will contain only real (float) coefficients. All corners must be at positive frequencies.
      COMPLEX ///< The filter *may* contain complex coefficients
    };

    /**
     * @class Corner
     *
     * @brief Specifies a point on a filter response curve. There is no
     * control of phase, just magnitude. 
     */
    struct Corner {
      /**
       * @brief constructor
       *
       * @param freq where the corner lands
       * @param gain the gain at that corner (in absolute, not dB)
       */
      Corner(float freq, float gain) : freq(freq), gain(gain) { }
      float freq;
      float gain;
    };

    /**
     * @class BadRealSpec
     *
     * @brief a subclass of std::runtime error thrown when a "real" valued filter spec is passed
     * a negative frequency.
     *
     * @param st a string explaining who was annoyed
     * @param freq the frequency at which it happened. 
     */
    class BadRealSpec : public std::runtime_error {
    public:
      BadRealSpec(const std::string & st, float freq) :
	std::runtime_error(SoDa::Format("FilterSpec::%1 specification for REAL valued filter contains a negative frequency %0. Not good.\n")
			 .addF(freq).addS(st).str()) { }
    };
    
    /**
     * @brief construct a filter specification
     *
     * @param sample_rate 
     * @param taps number of taps required
     * @param filter_type if REAL, then the filter shape (specified
     * by the added corners) must contain only positive
     * frequencies. The constructor will throw Filter::BadRealSpec
     * otherwise.
     */
    FilterSpec(float sample_rate, unsigned int taps,
	       FType filter_type = COMPLEX);



    /**
     * @brief Alternate constructor, for very simple band-pass filters
     * 
     * @param sample_rate in Hz -- as are all specified frequencies
     * @param low_cutoff lower 3dB point
     * @param high_cutoff upper 3dB point
     * @param skirt_width width of transition band
     * @param filter_type if REAL, then the filter shape (specified
     * by the added corners) must contain only positive
     * frequencies. The constructor will throw Filter::BadRealSpec
     * otherwise.
     * @param stop_band_attenuation_dB in dB
     */
    FilterSpec(float sample_rate, float low_cutoff, float high_cutoff, float skirt_width, 
	       FType filter_type = COMPLEX,
	       float stop_band_attenuation_dB = 60.0
	       );
    
    /**
     * copy constructor
     *
     */
    FilterSpec(const FilterSpec & f);

    /**
     * Set the starting gain for this filter. (Defaults to 1e-20)
     * @param gain the ideal amplitude at the lowest filter frequency (- sample_freq / 2) (absolute, not dB)
     * @return reference to this FilterSpec
     */
    FilterSpec & start(float gain);

    /**
     * Add a point in the transfer function. This specifies
     * a "corner" in the frequency response. 
     * 
     * @param freq the corner frequency
     * @param gain the ideal gain at the corner frequency (absolute, not dB)
     * @return reference to this FilterSpec
     */
    FilterSpec & add(float freq, float gain); 
      
    
    /**
     * @brief return a list of the corners in this filter
     * @return a list of the corners in this filter. 
     */
    const std::list<Corner> & getSpec();

    /**
     * @brief return the sample rate for which this filter will be created.
     *
     * @return the sample rate for which this filter will be created.
     */
    float getSampleRate() { return sample_rate; }

    /**
     * @brief build the frequency domain image of the prototype filter
     *
     * @param Hproto a buffer that will contain the DFT of the filter
     */
    void fillHproto(std::vector<std::complex<float>> & Hproto);

  protected:
    /**
     * @brief find the position in the prototype filter buffer
     *
     * @param freq the frequency who's bucket we are looking for
     * @return the bucket corresponding to freq
     */
    unsigned int indexHproto(float freq); 

  public:
    /**
     * @brief How many taps do we need to provide the requested 
     * shortest transition edge? 
     * 
     * Estimate and set the taps as appropriate. But the 
     * number of taps must be in the range min...max
     * 
     * @param min imum number of taps provided
     * @param max imum number of taps provided
     * @return number of taps chosen
     */
    unsigned int estimateTaps(unsigned int min, unsigned int max); 

    /**
     * @brief set the number of taps in the filter
     * @param new_tapcount yup...
     * @return the old tap count
     */
    unsigned int setTaps(unsigned int new_tapcount) { 
      taps = new_tapcount;
      return taps;
    }

    /**
     * @brief return the number of taps in the filter image
     *
     */
    unsigned int getTaps() {
      return taps; 
    }

    /**
     * @brief how wide is the filter?
     *
     * @return the frequency for the lowest and highest specified corner
     */
    std::pair<float, float> getFilterEdges();
    
  protected:
    bool sorted; ///< if true, the corner list is in order of frequency.
    /**
     * @brief sort the corner list from lowest frequency to highest
     */
    void sortSpec();
      
    unsigned int taps; 
    float sample_rate;
    FType filter_type;

    float stop_band_attenuation_dB; 
    
    std::list<Corner> spec; 
  };
}

