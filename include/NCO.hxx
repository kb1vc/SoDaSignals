#pragma once
/*
  Copyright (c) 2022 Matthew H. Reilly (kb1vc)
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

#include <complex>
#include <vector>

namespace SoDa {
  /**
   * 
   * @class NCO
   * Numerically controlled oscillator
   *
   * Frequency, sample rate, and all angles are set and reported in
   * double precision.
   *
   * Phase angle resolution has a lot to do with the cleanliness of
   * the oscillator output, even if the output type is complex<float>.
   */
  class NCO {
  public:
    /**
     * @brief Constructor
     *
     * @param sample_rate sample rate for the output stream
     * @param frequency frequency of the output stream
     */
    NCO(double sample_rate, double frequency);
    NCO();

    /**
     * @brief set the sample rate (change it)
     * @param sr new sample rate
     */
    void setSampleRate(double sr) { 
      sample_rate = sr; 
    }

    /**
     * @brief tune the NCO
     *
     * @param frequency make an output at this frequency
     */
    void setFreq(double frequency);

    /**
     * @brief adjust the current angle -- for phase alignment
     *
     * @param ang new angle in radians
     */
    void setAngle(double ang) { cur_angle = ang; }

    /**
     * @brief report the angle that will be used at the *next* step
     *
     * @return current phase angle in radians
     *
     */
    double getAngle() { return cur_angle; }

    /**
     * @brief report the amount of phase advance on each tick
     *
     * @return phase advance per tick (2 pi / freq)
     *
     */
    double getAngleIncr() { return ang_incr; }
    
    enum SumIt { 
      ADD, ///< get() method adds generated samples to the output vector
      SET  ///< get() method writes generated samples to the output vector
    };

    /**
     * @brief fill an length N vector with the next N steps of the oscillator
     *
     * @param out get will fill this vector (float)
     * @param sum if set to "ADD" the new steps will be *added* to the output, otherwise
     * the output value will be over-written
     *
     */
    void get(std::vector<std::complex<float>> & out, SumIt sum = SET);

    /**
     * @brief fill an length N vector with the next N steps of the oscillator
     *
     * @param out get will fill this vector (double)
     * @param sum if set to "ADD" the new steps will be *added* to the output, otherwise
     * the output value will be over-written
     *
     */
    void get(std::vector<std::complex<double>> & out, SumIt sum = SET);     

    /**
     * @class FreqOutOfBounds
     *
     * @brief the request would set the oscillator frequency above or below the nyquist limit for
     * this sample rate.
     */ 
    class FreqOutOfBounds : public std::runtime_error {
    public:
      /**
       * @brief Signal a violation of the nyquist limit
       *
       * @param fs the sample rate (Hz)
       * @param fr the requested frequency (Hz)
       */
      FreqOutOfBounds(double fs, double fr);
    };
    
  protected:
    double sample_rate; 
    double cur_angle;
    double ang_incr; 
  }; 
}

