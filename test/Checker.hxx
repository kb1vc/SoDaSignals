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
#include <vector>
#include <complex>
#include "../include/Periodogram.hxx"
#include "../include/NCO.hxx"

namespace SoDa {
  typedef std::vector<std::complex<float>> CVec; 
  class Checker {
  public:
    Checker(double sample_freq, uint32_t num_buckets = 2048);
    
    enum CheckRegion { STOP_BAND, PASS_BAND, TRANSITION_BAND };

    void setFreqAndReset(double freq, CheckRegion region); 

    void checkResponse(const CVec & data);
    
    void checkFinal();
    
    bool testPassed() { return test_passed && (sample_count > 32 * 1024); }

    uint32_t getNumBuckets();

    double getBucketWidth();

  protected:
    double bucketToFreq(uint32_t bucket);
    
    uint32_t num_buckets; 
    CheckRegion check_region; 
    double cur_freq;
    double sample_freq; 
    int pass_count;
    uint32_t sample_count; 
    SoDa::Periodogram * pdg_p;
    bool test_passed; 
    NCO ref_nco; 
  };
}
