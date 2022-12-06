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
      
    chk.setFreqAndReset(freq, reg); 

    osc.setFreq(freq); 

    for(int i = 0; i < num_passes; i++) {
      osc.get(in); 
      filt.apply(in, out); 
      chk.checkResponse(out); 
      if(i == 2) {
	// calculate the total power in and out
	float tot_in = 0.0;
	float tot_out = 0.0;
	for(auto v : in) {
	  auto av = std::abs(v);
	  tot_in += av * av; 
	}
	for(auto v : out) {
	  auto av = std::abs(v);
	  tot_out += av * av; 
	}
	std::cerr << SoDa::Format("RESP %0 %1 %2\n")
	  .addF(freq, 'e').addF(tot_in / in.size(), 'e')
	  .addF(tot_out / out.size(), 'e');
	  
      }
    }
    chk.checkFinal(); 
  
    if(!chk.testPassed()) {
      passed = false; 
      //      return false; 
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
