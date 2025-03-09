#pragma once

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

