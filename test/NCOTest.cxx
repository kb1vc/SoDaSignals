#include "../include/NCO.hxx"
#include <SoDa/Format.hxx>
#include <complex>
#include <vector>
#include <iostream>
#include <fstream>

int main() {
  SoDa::NCO nco(48e3, -500.0);

  std::vector<std::complex<float>> sig(1024);

  std::ofstream of("NCOtest.dat");
    
  for(int i = 0; i < 20; i++) {
    nco.get(sig); 
    for(auto v : sig) {
      of << SoDa::Format("%0 %1\n").addF(v.real(), 'e').addF(v.imag(), 'e');
    }
  }
  
  of.close();
  return 0;
}
