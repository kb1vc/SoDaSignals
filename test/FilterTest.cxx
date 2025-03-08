#include "../include/Filter.hxx"
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
int main() {

  float Fs = 48e3;
  float flp=-2.0e3;
  float fhp = 10e3;

  uint32_t buflen = 16 * 1024; 
  float skirt = 1000.0;
  float att = 50.0;
  uint32_t taps = int((Fs / skirt) * (att / 22.0));
  taps = taps | 1;
  SoDa::Filter filt(flp, fhp, 1e3, Fs, taps, buflen);

  std::cerr << "Taps = " << taps << "\n";
  
  struct amprec {
    amprec(float f, float p) { freq = f; power = p; }
    float freq; 
    float power; 
  };
  std::list<amprec> results; 
  
  for(float freq = -Fs / 2; freq < Fs / 2; freq += Fs * 0.001) {
    SoDa::NCO osc(Fs, freq); 
    std::vector<std::complex<float>> in(buflen), out(buflen); 
    osc.get(in); 
    filt.apply(in, out); 
    float mag = std::abs(out[buflen >> 1]);
    float pow = mag * mag;
    results.push_back(amprec(freq, pow)); 
  }
  
  bool is_ok = true; 
  // now look through the results...
  for(auto r : results) {
    auto logpow = 10 * std::log10(r.power);     
    if((r.freq < (flp - 3 * skirt)) || (r.freq > (fhp + 3 * skirt))) {
      // stop band. 
      if(logpow > -40) {
	std::cout << SoDa::Format("Out of band failure at freq %0 power level %1\n")
	  .addF(r.freq, 'e')
	  .addF(logpow, 'e');
	is_ok = false; 	
      }
    }
    else if((r.freq > flp + skirt * 0.25) && (r.freq < fhp - skirt * 0.25)) {
      // passband
      if((logpow < -1.0) || (logpow > 1.0)) {
	std::cout << SoDa::Format("In band failure at freq %0 power level %1\n")
	  .addF(r.freq, 'e')
	  .addF(logpow, 'e');
	is_ok = false; 		
      }
    }
    else {
      // transition band. 
      if(logpow > 0.0) {
	std::cout << SoDa::Format("Transition band failure at freq %0 power level %1\n")
	  .addF(r.freq, 'e')
	  .addF(logpow, 'e');
	is_ok = false; 		
      }
    }
  }
  
  if(is_ok) {
    std::cout << "PASSED\n";
  }
  else {
    std::cout << "FAILED\n";
  }
}
