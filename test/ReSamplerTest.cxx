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
    auto mag = std::abs(out[out.size() / 2]);
    auto pow = 20 * std::log10(mag);
    if((pow > 1.0) || (pow < -1.0)) {
      dumpVec("pdgout_weak_peak.dat", 
	      SoDa::Format("FAILED bad peak power: peak freq = %0 peak mag %1 peak_out %2 sf %3\n")
	      .addF(peak_freq, 'e').addF(mag, 'e').addF(peak_out, 'e').addF(sf, 'e').str(),
	      pdg_out, fs, freq);
      std::cout.flush();
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

  
  
  if(seg == OUT_OF_BAND) {
    for(int i = 0; i < pdg_out.size(); i++) {
      // are we below a respectable stop-band level for the
      // range?
      // is the level below -40db
      float vl; 
      if(pdg_out[i] < 1e-40) vl = -120; 
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

  uint32_t lx, ly; 
  resamp.getLxLy(lx, ly);
  std::cout << SoDa::Format("Resample %0 to %1 Lx = %2 Ly = %3  Filter len = %4 in size %6 out size %5\n")
    .addF(fs_in, 'e')
    .addF(fs_out, 'e')
    .addI(lx)
    .addI(ly)
    .addI(resamp.getFilterLength())
    .addI(resamp.getOutputBufferSize())
    .addI(resamp.getInputBufferSize());
  std::cout.flush();
  
  double fs_lim = std::max(fs_in, fs_out);
  SoDa::NCO nco(fs_in, fs_in * 0.125); // initial setting.
  SoDa::NCO ref_nco(fs_out, 0.0);
  auto out_buf_size = resamp.getOutputBufferSize();
  auto in_buf_size = resamp.getInputBufferSize();
  CVec in_buffer(in_buf_size);
  CVec out_buffer(out_buf_size); 
  CVec ref_out_buffer(out_buf_size); 
  
  uint32_t pdg_seg_len = 2048;
  SoDa::Periodogram pdg(pdg_seg_len); 
  
  double hfs_out = std::min(fs_in,fs_out) * 0.5;

  uint32_t num_sweeps = 4; 

  SoDa::Periodogram inpdg(pdg_seg_len);
  

  // are we inband? 
  BandSeg band_seg, last_band_seg; 
  last_band_seg = NONE; 


  // sweep the frequency. Do 1023 sample frequencies.
  double fd_increment = (fs_in / double(0.5 * pdg_seg_len));

  pdg.clear();
  int i = 0;
  int region = 0;
  for(double freq = -fs_in * 0.5; freq < fs_in * 0.5; freq += fd_increment) { // fd_increment) {
    bool first_phase_fail = false; // we haven't had a phase failure yet.

    i++;
    nco.setFreq(freq);
    

    pdg.clear();      
    
    auto fa = fabs(freq); 
    if(fa < 0.8 * hfs_out) {
      // in the passband.
      band_seg = PASS;
    }
    else if(fa < 1.1 * hfs_out) {
      band_seg = TRANSITION;
    }
    else {
      band_seg = OUT_OF_BAND;
    }

    
    last_band_seg = band_seg;
    
    if(band_seg == PASS) {
      ref_nco.setFreq(freq);
    }

    
    for(int i = 0; i < num_sweeps; i++) {
      // fill with the tone
      nco.get(in_buffer); 
      ref_nco.get(ref_out_buffer); 

      inpdg.accumulate(in_buffer);
      // resample
      resamp.apply(in_buffer, out_buffer);

      if((i == 2) && (band_seg == PASS)) {
	// sync the reference NCO
	// find the phase of the -10th element in the output buffer
	float ang = std::arg(out_buffer[out_buffer.size() >> 1]);
	float ref_ang = std::arg(ref_out_buffer[ref_out_buffer.size() >> 1]);
	// now align the reference NCO with the output stream.
	float nco_ang = ref_nco.getAngle();
	float ang_cor = ang - ref_ang; 
	fixAngle(ang_cor);
	ref_nco.setAngle(nco_ang + ang_cor);
	continue;
      }
      // accumulate the PDG
      // skip the first resample buffer. 
      // calculate the amplitude
      calcMag(in_buffer, out_buffer);

      pdg.accumulate(out_buffer);

      if(is_ok && (i > 2) && (band_seg == PASS)) {
	// check to see that the output buffer is smooth and
	// looks like the reference NCO.
	// don't worry about amplitude, just compare the angle;
	double ai = std::fabs(ref_nco.getAngleIncr());
	for(int j = 0; j < out_buffer.size(); j++) {
	  
	  auto ob_ang = std::arg(out_buffer[j]); 
	  auto ref_ang = std::arg(ref_out_buffer[j]); 
	  //fixAngle(ob_ang);
	  //fixAngle(ref_ang); 
	  float diff = ref_ang - ob_ang; 
	  fixAngle(diff);
	  if((std::fabs(diff) > ai) && (j > 20) && !first_phase_fail && (ai > 1e-5)) {
	    dump2CVec("badangle.dat", out_buffer, ref_out_buffer);
	    std::cout << "FAILED\n";
	    std::cout << SoDa::Format("Ooops... bad angle. freq %0 index %1 ref_ang = %2 ob angle %3  diff %4 bigger than %5 pass %6  ob (%7, %8)  ref (%9, %10)\n")
	      .addF(freq, 'e').addI(j)
	      .addF(ref_ang, 'e').addF(ob_ang, 'e')
	      .addF(diff, 'e').addF(ai, 'e')
	      .addI(i)
	      .addF(out_buffer[j].real(), 'e')
	      .addF(out_buffer[j].imag(), 'e')
	      .addF(ref_out_buffer[j].real(), 'e')
	      .addF(ref_out_buffer[j].imag(), 'e');
	    first_phase_fail = true;
	    std::vector<float> opdg;
	    pdg.get(opdg);
	    dumpVec("badanglePDG.dat", "whatever.", opdg, fs_in, freq);
	    is_ok = false; 
	    break; 
	  }
	}
      }
    }
    
    // check the pdg.
    if(!checkPDG(pdg, out_buffer, fs_out, freq, band_seg)) {
      std::cout << SoDa::Format("Bad result at freq %0 band seg %1\n").addF(freq, 'e').addI(band_seg);
      std::cout.flush();
      std::vector<float> ipdg; 
      inpdg.get(ipdg);
      dumpVec("inpdg.dat", "whatever.", ipdg, fs_in, freq);
      is_ok = false; 
    }
  }    
    //exit(-1);

  return is_ok;
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

  bool is_ok = true;

  typedef std::pair<double, double> dpair; 
  std::list<dpair> test_freqs = {
    dpair(48e3, 8e3),
    dpair(625e3, 48e3),
    dpair(120e3, 48e3),
    dpair(1.2e6, 48e3),
    dpair(1.25e6, 48e3),
    dpair(2.5e6, 48e3)
  };

  if(cmd.isPresent("fsin") && cmd.isPresent("fsout")) {
    std::cerr << SoDa::Format("Test %0 -> %1\n").addF(fs_in, 'e').addF(fs_out, 'e');
    is_ok = is_ok && testRatio(fs_in, fs_out);
  }
  else {
    for(auto fp : test_freqs) {
      std::cerr << SoDa::Format("Test %0 -> %1\n").addF(fp.first, 'e').addF(fp.second, 'e');
      is_ok = is_ok && testRatio(fp.first, fp.second);
      std::cerr << SoDa::Format("Test %1 -> %0\n").addF(fp.first, 'e').addF(fp.second, 'e');    
      is_ok = is_ok && testRatio(fp.second, fp.first);    
    }
  }

  
  if(is_ok) {
    std::cout << "\nPASSED\n";
  }
  else {
    std::cout << "\nFAILED\n";    
  }
}
