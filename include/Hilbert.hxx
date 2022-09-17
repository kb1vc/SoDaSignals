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


///
///  @file Hilbert.hxx
///
///  @brief This class creates a "hilbert transformer" to produce an
/// analytic (complex valued) signal from a real signal.
///
/// See Lyons (Understading Digital Signal Processing) pp 424... for a start
///  
///  @author M. H. Reilly (kb1vc)
///  @date   July 2022
///

#include <memory>
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include "FFT.hxx"
#include "OSFilter.hxx"
#include <stdexcept>

namespace SoDa {

  class Hilbert : public OSFilter{
  public:
    
    /// constructor
    /// Build the transformer for a particular buffer size.
    /// 
    /// @param taps number of taps in the impulse response. 
    /// @param in_buffer_size the input buffer will be this long

    Hilbert(unsigned int taps, unsigned int in_buffer_size); 

    ~Hilbert() {
      debug_of.close();
    }
    
    /// run the filter on a complex input stream
    /// @param in_buf the input buffer I/Q samples (real)
    /// @param out_buf the output buffer I/Q samples (complex)
    /// @return the length of the input buffer
    unsigned int apply(std::vector<float> & in_buf, 
		       std::vector<std::complex<float>> & out_buf);
    
    /**
     * @brief how big are the input/output vectors?
     *
     * @return in_buffer_size; 
     */
    unsigned int getSize() { return buffer_size; }

  protected:
    // the vector that holds our input real converted to complex
    std::vector<std::complex<float>> temp_in_buf;
    std::vector<std::complex<float>> temp_out_buf;    

    std::ofstream debug_of;
    
    unsigned int in_buffer_size; 
  };
}

