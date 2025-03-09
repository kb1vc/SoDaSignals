#pragma once
/*
 *
 *  @copyright  2025 Matt Reilly - kb1vc  All rights reserved.
 *
 *  @par BSD 2-Clause License
 *
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice, this
 *  list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <vector>
#include <complex>

namespace SoDa {

  /**
   * @brief calculate correlation between two vectors
   *
   * @param v0 first vector
   * @param v1 second vector
   *
   * @return the complex correlation (can be turned into mag and phase)
   * sum(v0 * conj(v1))
   */
  template<typename T>
  std::complex<T> correlate(const std::vector<std::complex<T>> & v0,
			    const std::vector<std::complex<T>> & v1) {
    std::complex<T> result(0, 0);

    auto vsize = (v0.size() > v1.size()) ? v1.size() : v0.size();    

    for(int i = 0; i < vsize; i++) {
      result += v0[i] * std::conj(v1[i]); 
    }

    result = result / T(vsize);

    return result;
  }

  /**
   * @brief calculate correlation between two vectors
   *
   * @param v0 first vector
   * @param v1 second vector
   *
   * @return the correlation magnitude
   */
  template<typename T>
  T correlate(const std::vector<T> & v0,
	      const std::vector<T> & v1) {
    T result = 0.0;

    auto vsize = (v0.size() > v1.size()) ? v1.size() : v0.size();

    for(int i = 0; i < vsize; i++) {
      result += v0[i] * v1[i];
    }

    T tvsize = vsize;
    result = result / tvsize;

    return result;
  
  }
  
}

