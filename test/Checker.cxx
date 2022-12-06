#include "Checker.hxx"
#include <iostream>
#include <fstream>
#include "../include/NCO.hxx"
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
namespace SoDa {
  typedef std::vector<std::complex<float>> CVec;
  typedef std::vector<float> FVec;   

  float fixAngle(float v) {
    while(v > M_PI) v = v - 2.0 * M_PI; 
    while(v < -M_PI) v = v + 2.0 * M_PI; 
    return v; 
  }
  
  void dump2CVec(std::ostream & of, const CVec & a, const CVec & b) {
    
    for(int i = 0; i < std::min(a.size(), b.size()); i++) {
      of << SoDa::Format("%0 %1 %2 %3 %4\n")
	.addI(i)
	.addF(a[i].real(), 'e')
	.addF(a[i].imag(), 'e')
	.addF(b[i].real(), 'e')
	.addF(b[i].imag(), 'e');
    }
  }

  void dump2CVec(const std::string & fname, const CVec & a, const CVec & b) {
    std::ofstream of(fname); 
    dump2CVec(of, a, b); 
    of.close();
  }
  
  
  void dumpFVec(std::ostream & of, const FVec & a) {
    for(auto v : a) {
      of << SoDa::Format("%0\n").addF(v, 'e');
    }
  }

  void dumpFVec(const std::string & fname, const FVec & a) {
    std::ofstream of(fname); 
    dumpFVec(of, a); 
    of.close();
  }
  
  Checker::Checker(double sample_freq, uint32_t num_buckets) : 
    sample_freq(sample_freq), num_buckets(num_buckets) {
    pass_count = 0; 
    ref_nco.setSampleRate(sample_freq); 
    pdg_p = new SoDa::Periodogram(num_buckets); 
    test_passed = true;
  }
  
  uint32_t Checker::getNumBuckets() { return num_buckets; }

  double Checker::getBucketWidth() { return sample_freq / double(num_buckets); }

  void Checker::setFreqAndReset(double freq, CheckRegion reg) {
    ref_nco.setFreq(freq); 
    pdg_p->clear();
    pass_count = 0; 
    sample_count = 0;
    test_passed = true; 
    check_region = reg;
    cur_freq = freq; 
  }

  void Checker::checkResponse(const CVec & data) {
    
    // cycle the NCO
    std::vector<std::complex<float>> osc(data.size());
    ref_nco.get(osc);
    
    sample_count += data.size(); 
    pass_count++;
    if(pass_count < 2) return;

    // add to the periodogram 
    pdg_p->accumulate(data); 


    if((pass_count == 2) && (check_region == PASS_BAND)){
      // do the phase alignment, but only if we are in the passband. 
      float ang = std::arg(data[data.size() - 10]);
      float ref_ang = std::arg(osc[data.size() - 10]);
      // align the reference oscilator
      float nco_ang = ref_nco.getAngle();// + ref_nco.getAngleIncr();;
      float ang_cor = ang - ref_ang;
      ref_nco.setAngle(fixAngle(nco_ang + ang_cor));
      return;
    }

    auto ai = ref_nco.getAngleIncr();

    // do the phase checks
    if(check_region == PASS_BAND) {
      for(int i = 0; i < data.size(); i++) {
	float ang_in = std::arg(data[i]);
	float ang_ref = std::arg(osc[i]);
	float diff = std::abs(ang_in - ang_ref); 
	fixAngle(diff); 
	float target = std::abs(ai) * 0.1;
	bool is_ok = (diff < ai) || 
	  ((std::abs(diff) - (2 * M_PI)) < ai);
	if(!is_ok) {
	  if(test_passed) {
	    // we're the cause for the first failure. 
	    std::cerr << SoDa::Format("Bad inband phase at freq %0, index %1 diff %2 target %3 ang_in %4 ang_ref %5\n")
	      .addF(cur_freq, 'e')
	      .addI(i)
	      .addF(diff, 'e')
	      .addF(std::abs(ai) * 0.1)
	      .addF(ang_in)
	      .addF(ang_ref); 
	    dump2CVec("PhaseFailure.dat", data, osc);
	  }
	  test_passed = false; 
	}
      }
    }

    return; 
  }

  double Checker::bucketToFreq(uint32_t bucket) {
    int psize = int(pdg_p->getSize());
    int bkt = int(bucket);
    
    double bucket_width = sample_freq / double(psize);
    double ret =  double(bkt - (psize / 2)) * bucket_width; 
    std::cerr << SoDa::Format("bucket %0 bw = %1 ret %2 pdg size %3\n")
      .addI(bkt)
      .addF(bucket_width, 'e')
      .addF(ret, 'e')
      .addI(psize);
    return ret; 
  }

  void Checker::checkFinal() {
    // now we do all the check against the periodogram

    // frequency in the passband? is it correct? is the phase good?
    // find the peak
    
    std::vector<float> pdg_out;
    pdg_p->get(pdg_out); 
    auto pdg_sf = pdg_p->getScaleFactor();
    for(auto & v : pdg_out) {
      v = v * pdg_sf;
    }

    double bucket_size = sample_freq / pdg_out.size();
    
    // in the stop band?  is everybody below the threshold?
    if(check_region == STOP_BAND) {
      // is everybody below -40 dB?
      // that would be 1e-20 in abs magnitude
      int bucket = 0; 
      for(auto v : pdg_out) {
	if(std::log10(std::abs(v)) > -2.0) {
	  if(test_passed) {
	    // we're the cause for the first failure. 
	    std::cerr << SoDa::Format("Bad out-of-band power spike! at bucket %0 (freq %1) power %2 (%3 dB)\n")
	      .addI(bucket)	      
	      .addF(bucketToFreq(bucket), 'e')
	      .addF(std::abs(v), 'e')
	      .addF(20.0 * std::log10(std::abs(v)), 'e');

	    dumpFVec("MagnitudeFailure.dat", pdg_out);
	  }
	  test_passed = false; 
	}
	bucket++; 
      }
      return;
    }

    // otherwise we're in the passband or transition band.
    // find the peak and make sure it is us. 
    int bucket = 0; 
    int peak_bucket = 0;
    float peak_amp = 0.0;
    for(auto v : pdg_out) {
      
      if(std::abs(v) > peak_amp) {
	peak_amp = std::abs(v); 
	peak_bucket = bucket; 
      }
      bucket++; 
    }
   using the periodogram to test for magnitude/total power is a bad idea.
						 
    std::cerr << SoDa::Format("PB %0 %1\n").addF(cur_freq, 'e').addF(20 * std::log10(peak_amp));

    // is the peak near us?
    int correct_bucket = int((cur_freq + (sample_freq / 2))  / bucket_size + 0.5); 
    if(std::abs(peak_bucket - correct_bucket) > 1) {
      if(test_passed) {
	std::cerr << SoDa::Format("Bad peak frequency found at bucket %0 freq %1  should be %2 freq %3\n")
	  .addI(peak_bucket)
	  .addF(bucketToFreq(peak_bucket), 'e')
	  .addI(correct_bucket)
	  .addF(cur_freq, 'e');
	dumpFVec("FrequencyFailure.dat", pdg_out);
	test_passed = false; 
      }
    }

    int cor_peak_bucket = peak_bucket - (pdg_out.size() / 2);
    // in the transition band? Are we below 0dB? 
    if(check_region == PASS_BAND) {
      if((peak_amp < 0.8) || (peak_amp > 1.2)) { // more than 1dB ripple
	if(test_passed) {
	  // only print on the first failure. 
	  std::cerr << SoDa::Format("Bad pass band response found  bucket %0 (cor %3) freq %1  mag was %2 dB\n")
	    .addI(peak_bucket)
	    .addF(bucketToFreq(cor_peak_bucket), 'e')
	    .addF(std::log10(std::abs(peak_amp)) * 20.0, 'e')
	    .addI(cor_peak_bucket);	  
	  dumpFVec("PassBandFailure.dat", pdg_out);
	  test_passed = false; 
	}
      }
    }
    else {
      // must be transition. 
      if(peak_amp > 1.0) {
	if(test_passed) {
	  std::cerr << SoDa::Format("Bad transition band response found  bucket %0 freq %1  mag was %2 dB\n")
	    .addI(peak_bucket)
	    .addF(bucketToFreq(peak_bucket), 'e')
	    .addF(std::log10(std::abs(peak_amp) * 20.0), 'e');
	  dumpFVec("TransitionBandFailure.dat", pdg_out);
	  test_passed = false; 
	}
      }
    }
  }
}

