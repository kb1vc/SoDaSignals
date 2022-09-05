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
  float Fs = 125e3;
  float flp = -20e3;
  float fhp = 40e3;

  uint32_t buflen = 16 * 1024; 
  float skirt = 2e3;
  float att = 50.0;
  uint32_t taps = int((Fs / skirt) * (att / 22.0));
  taps = taps | 1;
  taps = 53;
  Fs = 48e3;
  flp=-2.0e3;
  fhp = 10e3;
  SoDa::Filter filt(flp, fhp, 1e3, Fs, taps, buflen);

  std::cerr << "Taps = " << taps << "\n";
  
  for(float freq = -Fs / 2; freq < Fs / 2; freq += Fs * 0.001) {
    SoDa::NCO osc(Fs, freq); 
    std::vector<std::complex<float>> in(buflen), out(buflen); 
    osc.get(in); 
    filt.apply(in, out); 
    std::cout << SoDa::Format("%0 %1 %2 %3\n").addF(freq, 'e').addF(std::abs(in[buflen >> 1]))
      .addF(std::abs(out[buflen >> 1]))
      .addF(fixAngle(std::arg(in[buflen >> 1]) - std::arg(out[buflen >> 1])), 'e');

  }
}
