#include <iostream>
#include <fstream>
#include <chrono>

#include "../include/Filter.hxx"
#include "../include/OSFilter.hxx"
#include "../include/NCO.hxx"
#include "Checker.hxx"
#include <SoDa/Format.hxx>
#include <cmath>

float fixAngle(float a) {
  while(a > M_PI) a = a - 2.0 * M_PI;
  while(a < -M_PI) a = a + 2.0 * M_PI; 
  return a; 
}

SoDa::Checker::CheckRegion getRegion(double freq, double flo, double fhi, double skirt) {
  SoDa::Checker::CheckRegion reg; 
  if((freq > flo + 0.95 * skirt) && (freq < (fhi - 0.95 * skirt))) {
    reg = SoDa::Checker::PASS_BAND; 
    std::cerr << SoDa::Format("Test freq %0 flo %1 fhi %2 skirt %3\n")
      .addF(freq, 'e').addF(flo, 'e').addF(fhi, 'e').addF(skirt, 'e');
  }
  else if((freq < (flo - 1.05 * skirt)) || (freq > (fhi + 1.05 * skirt))) {
    reg = SoDa::Checker::STOP_BAND;       
    std::cerr << SoDa::Format("Test freq %0 flo %1 fhi %2 skirt %3\n")
      .addF(freq, 'e').addF(flo, 'e').addF(fhi, 'e').addF(skirt, 'e');
  }
  else {
    reg = SoDa::Checker::TRANSITION_BAND;             
  }

  return reg; 
}

bool test() {
  double Fs = 48.0e3;
  double flo = -2.0e3;
  double fhi = 10.0e3;
  uint32_t buflen = 8304;
  float skirt = 1.0e3;

  SoDa::OSFilter filt(flo, fhi, skirt, Fs, buflen);
  
  auto edges = filt.getFilterEdges(); 
  std::cerr << SoDa::Format("Filter taps %0 fft size %1 edges %2 %3\n")
    .addI(filt.getTaps())
    .addI(filt.getInternalSize())
    .addF(edges.first, 'e')
    .addF(edges.second, 'e');
  
  int num_passes = 6;

  SoDa::Checker chk(Fs, filt.getTaps(), 1.0, 35.0, 0.01, buflen, 256 * 4); 

  bool passed = true; 

  for(uint32_t i = 0; i < chk.getNumFreqSteps(); i++) {
    chk.checkResponse(i, [flo, fhi, skirt](double f) { return getRegion(f, flo, fhi, skirt); },
		      [& filt](std::vector<std::complex<float>> & in,
			       std::vector<std::complex<float>> & out) { filt.apply(in, out); });
    
    if(!chk.testPassed()) {
      passed = false; 
    }
  }
  
  return passed;
}

int main() {
  if(test()) {
    std::cout << "PASSED\n";
  }
  else {
    std::cout << "FAILED\n";
  }
}
