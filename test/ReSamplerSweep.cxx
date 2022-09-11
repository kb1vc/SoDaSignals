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

#include "../include/ReSampler.hxx"
#include "../include/Periodogram.hxx"
#include "../include/NCO.hxx"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <SoDa/Options.hxx>
#include <SoDa/Format.hxx>

typedef std::vector<std::complex<float>> CVec;
typedef std::vector<float> FVec; 

// the test strategy:
//
// Put various frequencies into the resampler across the input band.
// then calculate a periodogram and look for a single spike at the
// appropriate frequency (if it is inband for the output) and
// a null if it is out of band in the output.
//
// This borrows some code from the Periodogram test... 

void fixAngle(float & a) {
  while(a > M_PI) a = a - 2.0 * M_PI;
  while(a < -M_PI) a = a + 2.0 * M_PI;
}

void dump2CVec(std::ostream & of, CVec a, CVec b) {

  for(int i = 0; i < std::min(a.size(), b.size()); i++) {
    of << SoDa::Format("%0 %1 %2 %3 %4\n")
      .addI(i)
      .addF(a[i].real(), 'e')
      .addF(a[i].imag(), 'e')
      .addF(b[i].real(), 'e')
      .addF(b[i].imag(), 'e');
  }
}
void dump2CVec(const std::string & fn, CVec a, CVec b) {
  std::ofstream of(fn);
  dump2CVec(of, a, b); 
  of.close();
}

void dumpVec(std::ostream & os, const std::string & msg, FVec v, double Fs, double Fr) {
  double freq_inc = Fs / double(v.size());
  std::cout << SoDa::Format("%0 Fs = %1 Fr = %2\n")
    .addS(msg)
    .addF(Fs, 'e')
    .addF(Fr, 'e');

  for(int i = 0; i < v.size(); i++) {
    os << SoDa::Format("%0 %1\n").addF(double(i) * freq_inc - Fs / 2, 'e').addF(v[i]);
  }
}
void dumpVec(const std::string fname, const std::string & msg, FVec v, double Fs, double Fr) {
  std::ofstream of(fname); 
  dumpVec(of, msg, v, Fs, Fr); 
  of.close();
}

void calcMag(CVec & in_buffer, CVec & out_buffer) {
  float imax = 0.0; 
  for(auto v : in_buffer) {
    auto mag = std::abs(v); 
    imax = std::max(mag, imax); 
  }
  float omax = 0.0;
  for(auto v : out_buffer) {
    auto mag = std::abs(v); 
    omax = std::max(mag, omax); 
  }
  // std::cout << SoDa::Format("Input mag %0 (%1 dB)  Output mag %2 (%3 dB)\n")
  //   .addF(imax, 'e')
  //   .addF(20.0 * std::log10(imax), 'e')
  //   .addF(omax, 'e')
  //   .addF(20.0 * std::log10(omax), 'e');
}

enum BandSeg { NONE, PASS, TRANSITION, OUT_OF_BAND }; 

bool checkPDG(SoDa::Periodogram & pdg, 
	      std::vector<std::complex<float>> & out, 
	      double fs,
	      double freq, 
	      BandSeg seg) {
  // get the periodogram
  std::vector<float> pdg_out; 
  pdg.get(pdg_out); 

  auto sf = pdg.getScaleFactor();

  // where is the peak index?
  int peak_idx = 0;
  float peak = 0.0;
  int idx = 0;
  for(auto & v : pdg_out) {
    v = v * sf; 
    if(v > peak) {
      peak = v; 
      peak_idx = idx; 
    }
    idx++; 
  }
  
  // is the peak input frequency in the peak bucket? 
  double bucket_size = (fs / double(pdg.getSize()));
  double peak_freq = bucket_size * peak_idx - 0.5 * fs;

  float peak_out = std::abs(pdg_out[peak_idx]);
  if(seg == PASS) {
    auto mag = std::abs(peak_out * sf);
    if((mag > 1.2) || (mag < 0.9)) {
      dumpVec("pdgout_weak_peak.dat", 
	      SoDa::Format("FAILED bad peak power: peak freq = %0 peak mag %1 peak_out %2 sf %3\n")
	      .addF(peak_freq, 'e').addF(mag, 'e').addF(peak_out, 'e').addF(sf, 'e').str(),
	      pdg_out, fs, freq);
      std::cout.flush();
      exit(-1);
      return false; 
    }
  }
  
  if((seg != OUT_OF_BAND) && 
     (fabs(freq - peak_freq) > bucket_size) &&
     (20.0 * std::log10(std::fabs(pdg_out[peak_idx])) > ((seg == OUT_OF_BAND) ? -40 : -2))) {
    dumpVec("pdgout_bad_peak.dat", SoDa::Format("Peak freq = %0 freq = %1 bucket_size = %2 mag = %3 dB\n")
	    .addF(peak_freq,'e').addF(freq, 'e').addF(bucket_size, 'e')
	    .addF(20 * std::log10(std::abs(pdg_out[peak_idx])), 'e').str(), 
	    pdg_out, fs, freq);
    return false; 
  }
  
  bool is_ok = true; 
  int low_rs_cutoff_idx = peak_idx - 3;
  int high_rs_cutoff_idx = peak_idx + 3;

  auto log_peak = std::log10(peak);

  
  
  for(int i = 0; i < pdg_out.size(); i++) {
    // are we below a respectable stop-band level for the
    // range?
    if((i < low_rs_cutoff_idx) || (i > high_rs_cutoff_idx)) {
      // is the level below -40db
      float vl; 
      if(pdg_out[i] < 1e-20) vl = -120; 
      else  vl = 20 * std::log10(pdg_out[i]);
      if(vl > -40) {
	is_ok = false;
	dumpVec("pdg_out_of_band.dat", SoDa::Format(R"(Out of band response freq = %0 freq = %1 bucket_size = %2
   i = %3  low_rs_cutoff_idx %4 high_rs_cutof_idx %5  pow %6
)")
		.addF(freq,'e').addF(peak_freq, 'e').addF(bucket_size, 'e')
		.addI(i).addI(low_rs_cutoff_idx).addI(high_rs_cutoff_idx)
		.addF(vl, 'e')
		.str(), 
		pdg_out, fs, freq);
	return false; 
      }
	
    }
  }
  
  return is_ok; 
}

bool testRatio(double fs_in, double fs_out) {
  bool is_ok = true;
  
  // first create the resampler sized for 50 mS blocks.
  SoDa::ReSampler resamp(fs_in, fs_out, 0.05); 

  auto ncosr = std::max(fs_out, fs_in);
  SoDa::NCO nco(ncosr, ncosr * 0.125); // initial setting.  
  double fs_lim = std::max(fs_in, fs_out);
  auto out_buf_size = resamp.getOutputBufferSize();
  auto in_buf_size = resamp.getInputBufferSize();
  CVec in_buffer(in_buf_size);
  CVec out_buffer(out_buf_size); 
  
  uint32_t pdg_seg_len = 1024; 
  SoDa::Periodogram pdg(pdg_seg_len); 
  

  double hfs_out = std::min(fs_in,fs_out) * 0.5;

  uint32_t num_sweeps = 4; 

  SoDa::Periodogram inpdg(pdg_seg_len);
  

  // are we inband? 
  BandSeg band_seg, last_band_seg; 
  last_band_seg = NONE; 


  // sweep the frequency. Do 1023 sample frequencies.
  double fd_increment = fs_in / double(1253);

  pdg.clear();
  int i = 0;
  // sweep the output passband from -fs_out/2 to fs_out/2
  double freq_incr = fs_out / double(pdg_seg_len);
  bool collect = false; 
  for(double freq = -fs_in * 0.5; freq < fs_in * 0.5; freq += fs_out / double(pdg_seg_len)) { 

    if ((freq > 0.0) && (freq < freq_incr * 1.1)) {
      if(collect == false) {
	collect = true; 
	pdg.clear();
      }
    }
    else {
      if(collect) {
	std::vector<float> pdg_out;
	auto sf = pdg.getScaleFactor();
	pdg.get(pdg_out);
	for(auto & v : pdg_out) {
	  v = v * sf; 
	}
	dumpVec("LowFreqPDG.dat", "whatever", pdg_out, fs_out, freq - freq_incr);
	collect = false; 
      }
    }

    nco.setFreq(freq);
    
    float total_power = 0.0; 
    for(int j = 0; j < num_sweeps; j++) {
      // fetch a new block
      nco.get(in_buffer); 
      // resample
      resamp.apply(in_buffer, out_buffer); 

      if(j > 0) {
	// accumulate total power in output stream
	for(int k = 0; k < out_buffer.size(); k++) {
	  float mag = std::abs(out_buffer[k]);
	  total_power += mag * mag; 
	}
	if(collect) {
	  pdg.accumulate(out_buffer); 
	}
	float  mag = std::abs(out_buffer[out_buffer.size() / 2]);
	total_power = mag*mag; 
      }
    }
    auto resamp_rat = (fs_out / fs_in);
    auto resamp_sf = resamp_rat * resamp_rat;
    
    // now print the average power vs. frequency
    std::cout << SoDa::Format("%0 %1\n")
      .addF(freq, 'e')
      //      .addF(total_power / float((num_sweeps - 1)), 'e'); 
      .addF(total_power, 'e');
  }
  return true; 
}


int main(int argc, char * argv[]) {
  int sample_rate = 1000; 
  // a bandpass filter from

  double fs_in;
  double fs_out; 
  
  SoDa::Options cmd;
  cmd.add(&fs_in, "fsin", 'i', 48e3, "interpolate")
    .add(&fs_out, "fsout", 'd', 8e3, "decimate")
    .addInfo("Test rational resampler with sweeping input.\n");

  if(!cmd.parse(argc, argv)) {
    std::cout << "Bad command line\nFAILED\n";
    exit(-1);
  }

  if(cmd.isPresent("fsin") && cmd.isPresent("fsout")) {
    std::cerr << SoDa::Format("Test %0 -> %1\n").addF(fs_in, 'e').addF(fs_out, 'e');
    testRatio(fs_in, fs_out);
  }
}
