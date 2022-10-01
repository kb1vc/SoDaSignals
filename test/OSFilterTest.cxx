#include "../include/Filter.hxx"
#include "../include/OSFilter.hxx"
#include "../include/NCO.hxx"
#include "Checker.hxx"
#include <SoDa/Format.hxx>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>

float fixAngle(float a) {
  while(a > M_PI) a = a - 2.0 * M_PI;
  while(a < -M_PI) a = a + 2.0 * M_PI; 
  return a; 
}

int orig() {
  float Fs = 125e3;
  float flp = -20e3;
  float fhp = 40e3;

  uint32_t buflen = 8304; 
  float skirt = 2e3;
  float att = 50.0;
  uint32_t taps = int((Fs / skirt) * (att / 22.0));
  taps = taps | 1;
  taps = 337;
  Fs = 48e3;
  flp=-2.0e3;
  fhp = 10e3;

  SoDa::OSFilter filt(flp, fhp, 1e3, Fs, buflen);

  std::cerr << "Taps = " << filt.getTaps() << "\n";

  int num_passes = 4; 

  for(float freq = -Fs / 2; freq < Fs / 2; freq += Fs * 0.01) {
    SoDa::NCO osc(Fs, freq); 
    std::vector<std::complex<float>> in(buflen), out(buflen); 
    for(int i = 0; i < num_passes; i++) {
      osc.get(in); 
      filt.apply(in, out);

      std::cout << SoDa::Format("%0 %1 %2 %3\n").addF(freq, 'e').addF(std::abs(in[buflen >> 1]))
	.addF(std::abs(out[buflen >> 1]))
	.addF(fixAngle(std::arg(in[buflen >> 1]) - std::arg(out[buflen >> 1])), 'e');
    }
    break; 
  }
  return 0; 
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

  SoDa::Checker chk(Fs); 

  bool passed = true; 

  std::vector<std::complex<float>> in(buflen), out(buflen); 
  SoDa::NCO osc(Fs, Fs / 4);
  
  auto bw = chk.getBucketWidth();

  for(double freq = -Fs/2 + (Fs/1024); freq < Fs/2; freq += Fs / (1.5 * buflen)) {
    SoDa::Checker::CheckRegion reg; 
    if((freq > flo + skirt) && (freq < (fhi - skirt))) {
      reg = SoDa::Checker::PASS_BAND; 
    }
    else if((freq < (flo - skirt)) || (freq > (fhi + skirt))) {
      reg = SoDa::Checker::STOP_BAND;       
    }
    else {
      reg = SoDa::Checker::TRANSITION_BAND;             
    }
      
    chk.setFreqAndReset(freq, reg); 

    osc.setFreq(freq); 
  
    for(int i = 0; i < num_passes; i++) {
      osc.get(in); 
      filt.apply(in, out); 
      chk.checkResponse(out); 
    }
    chk.checkFinal(); 
  
    if(!chk.testPassed()) {
      passed = false; 
      return false; 
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
