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

    void setFreq(double frequency);

    void setAngle(double ang) { cur_angle = ang; }
    double getAngle() { return cur_angle; }
    double getAngleIncr() { return ang_incr; }
    
    enum SumIt { ADD, SET };
    
    void get(std::vector<std::complex<float>> & out, SumIt sum = SET);

    void get(std::vector<std::complex<double>> & out, SumIt sum = SET);     

    class FreqOutOfBounds : public std::runtime_error {
    public:
      FreqOutOfBounds(double fs, double fr);
    };
    
  protected:
    double sample_rate; 
    double cur_angle;
    double ang_incr; 
  }; 
}

