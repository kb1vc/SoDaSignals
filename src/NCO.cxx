#include "../include/NCO.hxx"
#include <SoDa/Format.hxx>
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

namespace SoDa {
  NCO::NCO(double _sample_rate, double frequency)  {
    sample_rate = _sample_rate; 
    setFreq(frequency); 
    cur_angle = 0.0; 
  }
  
  void NCO::setFreq(double frequency) {
    if(fabs(frequency) > 0.5 * sample_rate) {
      throw FreqOutOfBounds(sample_rate, frequency); 
    }
    ang_incr = M_PI * 2.0 * frequency / sample_rate; 
  }
  
  template<typename T> 
  void genGet(std::vector<std::complex<T>> & out, 
	      double & angle, double ang_incr, 
	      NCO::SumIt sum) {
    for(auto & v : out) {
      if(sum == NCO::ADD) {
	v += std::complex<T>(cos(angle), sin(angle));
      }
      else {
	v = std::complex<T>(cos(angle), sin(angle));	
      }
      angle += ang_incr; 
      if(angle > M_PI) angle = angle - 2.0 * M_PI;
      if(angle < -M_PI) angle = angle + 2.0 * M_PI; 
    }
  }
  
  void NCO::get(std::vector<std::complex<float>> & out, SumIt sum) {
    genGet<float>(out, cur_angle, ang_incr, sum);
  }

  void NCO::get(std::vector<std::complex<double>> & out, SumIt sum) {
    genGet<double>(out, cur_angle, ang_incr, sum);
  }

  NCO::FreqOutOfBounds::FreqOutOfBounds(double fs, double fr) : 
    std::runtime_error(SoDa::Format("SoDa::NCO::setFreq Frequency is out-of bounds. Sample Freq %0 Requested Freq %1\n").addF(fs, 'e').addF(fr, 'e').str()) {
  }
} 
