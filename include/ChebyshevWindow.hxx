#pragma once
#include <vector>

namespace SoDa {
  void ChebyshevWindow(std::vector<float> & win, size_t N, float atten);
  void ChebyshevWindow(std::vector<double> & win, size_t N, float atten);   
}
