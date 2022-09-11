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

#include "FilterSpec.hxx"
#include <cmath>

namespace SoDa {
  FilterSpec::FilterSpec(float sample_rate, unsigned int taps, FType filter_type) :
    sorted(false), filter_type(filter_type), taps(taps), sample_rate(sample_rate)
  {
  }

  FilterSpec::FilterSpec(float sample_rate, unsigned int taps, float low_cutoff, float high_cutoff, 
			 float skirt_width, 
			 FType filter_type) :
    sorted(false), filter_type(filter_type), taps(taps), sample_rate(sample_rate)
  {
    int lo_idx = indexHproto(low_cutoff); 
    int hi_idx = indexHproto(high_cutoff);
    // now backwards figure the frequencies to add...
    float hz_per_bucket = sample_rate / float(taps); 
    float low_stop = low_cutoff -  1.0001 * hz_per_bucket;
    float high_stop = high_cutoff + 1.0001 * hz_per_bucket;
    float low_pass = low_cutoff + 1.0001 * hz_per_bucket; 
    float high_pass = high_cutoff - 1.0001 * hz_per_bucket;
    
    add(low_stop, 0.0);
    add(low_cutoff, 1.0);
    add(high_cutoff, 1.0);
    add(high_stop, 0.0);
  }
  
  unsigned int FilterSpec::estimateTaps(unsigned int min_taps, unsigned int max_taps) {
    // find the narrowest transition region
    if(!sorted) sortSpec();

    // now crawl...
    bool looking = true;
    float min_interval = 10e9;
    float last_gain;
    float last_corner;
    last_gain = spec.front().gain;
    last_corner = spec.front().freq;
    for(const auto & v : spec) {
      if(v.gain != last_gain) {
	float inter = v.freq - last_corner;
	if(inter < min_interval) min_interval = inter;
	last_gain = v.gain;
	last_corner = v.freq;
      }
    }
    
    // now use fred harris's rule for the number of filter taps?
    // assume 40 dB stop band attenuation
    float ftaps = sample_rate * 40 / (min_interval * 22);

    unsigned int comp_taps = 2 * int(ftaps / 2) + 1; // make sure it is odd

    taps = std::max(comp_taps, min_taps);
    taps = std::min(taps, max_taps);

    return taps;
  }

  void FilterSpec::fillHproto(std::vector<std::complex<float>> & Hproto) {
    if(!sorted) sortSpec();

    Hproto.resize(taps);
    
    std::list<std::pair<Corner,Corner>> edges;
    float start_freq = (filter_type == REAL) ? 0.0 : - (sample_rate / 2);
    float end_freq = sample_rate/2;
    Corner last(start_freq, -200.0);

    for(int i = 0; i < taps; i++) Hproto.at(i) = std::complex<float>(0.0,0.0);
    for(auto v : spec) {
      auto idx = indexHproto(v.freq); 
      for(int i = idx; i < taps; i++) {
	Hproto.at(i) = std::complex<float>(v.gain,0.0);
      }
    }
  }

  FilterSpec & FilterSpec::start(float amp) {
    sorted = false;
    float start_freq = (filter_type == REAL) ? 0.0 : - (sample_rate / 2);

    if(amp < -200) amp = -200;
    
    spec.push_back(Corner(start_freq, amp));
    
    return *this;
  }

  FilterSpec & FilterSpec::add(float freq, float amp) {
    sorted = false;
    if((filter_type == REAL) && (freq < 0)) {
      throw BadRealSpec("add", freq);
    }
    // no zero amplitude buckets. (We need to be able to divide and take logs.)
    if(amp < -200) amp = -200;
    
    spec.push_back(Corner(freq, amp));    
    
    return *this;
  }
      
  const std::list<FilterSpec::Corner> FilterSpec::getSpec() {
    if(!sorted) sortSpec();
    return spec; 
  }

  std::pair<float, float> FilterSpec::getFilterEdges() {
    std::pair<float, float> ret;
    ret.first = (filter_type == REAL) ? 0.0 : -(sample_rate / 2);
    ret.second = sample_rate / 2;
    bool looking_low;
    float last_freq; 
    for(const auto v : spec) {
      if(looking_low) {
	if(v.gain > (1 - 0.01)) {
	  ret.first = v.freq;
	  looking_low = false; 
	}
      }
      else {
	if(v.gain < 0.1) {
	  ret.second = last_freq;
	  break; 
	}
      }
      last_freq = v.freq; 
    }
    return ret; 
  }
  
  void FilterSpec::sortSpec() {
    spec.sort([](const Corner & a, const Corner & b) { return a.freq < b.freq; });
    sorted = true; 
  }

  unsigned int FilterSpec::indexHproto(float freq) {
    unsigned int ret;
    float hsamprate = sample_rate / 2;
    float hz_per_bucket = sample_rate / float(taps);
    float norm_freq = freq / hz_per_bucket;
    int bucket = int(norm_freq + 0.5);
    ret = ((taps - 1) / 2) + bucket;
    if(bucket >= 0) ret++;
    
    if(ret >= taps) ret = taps - 1;
    
    return ret;
  }
}
