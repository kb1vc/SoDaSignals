#include "Filter.hxx"
#include "ChebyshevWindow.hxx"
#include <stdexcept>
#include <sstream>
#include <iostream>

void SoDa::Filter::apply(std::vector<std::complex<float>> & out, std::vector<std::complex<float>> & in) {
  if(H.size() != in.size()) makeImage(in.size());
  
  // take the FFT of the input
  dft_p->fft(out, in);
  // now multiply by the filter image
  // and the gain correction
  float gain_corr = 1.0 / ((float) in.size());
  for(int i = 0; i < in.size(); i++) {
    temp[i] = out[i] * H[i] * gain_corr; 
  }

  // now the inverse fft
  dft_p->ifft(out, temp); 
}

void SoDa::Filter::applyCont(std::vector<std::complex<float>> & out, std::vector<std::complex<float>> & in) {
  int i, j;
  debug_count++;
  
  if(last_invec_size != in.size()) {
    setupNewOverlap(in.size()); 
  }

  // put the new samples at the end of the saved buffer.
  for(i = 0; i < in.size(); i++) {
    saved[i + save_length] = in[i];
  }
  
  // take the FFT of the saved buffer
  dft_p->fft(saved_dft, saved);

  // multiply by the image
  for(i = 0; i < saved_dft.size(); i++) {
    saved_dft[i] = H[i] * saved_dft[i];
  }
  
  // take the inverse fft
  dft_p->ifft(temp, saved_dft); 

  // copy the result to the output
  float scale = 1.0 / ((float) temp.size());
  for(i = 0; i < out.size(); i++) {
    out[i] = temp[i] * scale; 
  }

  // save the last section of the input
  for(i = 0, j = (in.size() - save_length); i < save_length; i++, j++) {
    saved[i] = in[j]; 
  }

  
}


// normalize a frequency (in Hz) to the range -pi, pi
float SoDa::Filter::normalize(float freq, float sample_rate) {
  float res = freq / sample_rate;

  res = M_PI * res;

  if(res > M_PI) res = res - 2.0 * M_PI;
  if(res < -M_PI) res = res + 2.0 * M_PI; 
  return res; 
}

void SoDa::Filter::makePrototype(std::vector<std::complex<float>> & PH, float sample_rate, float lo, float hi, 
				 float transition_width) {
  float lo_stop = lo - transition_width;
  float lo_pass = lo;
  float hi_pass = hi;
  float hi_stop = hi + transition_width;

  int num_taps = PH.size();
  float bucket_width = sample_rate / ((float) num_taps);

  for(int i = 0; i < num_taps; i++) {
    int f_idx = i - (num_taps / 2);
    float freq = bucket_width * ((float) f_idx); 
    float v; 
    if((freq <= lo_stop) || (freq >= hi_stop)) {
      PH[i] = 0.0;
    }
    else if((freq >= lo_pass) && (freq <= hi_pass)) {
      PH[i] = 1.0;
    }
    else if((freq > lo_stop) && (freq < lo_pass)) {
      v = (freq - lo_stop) / (lo_pass - lo_stop);
      PH[i] = v; 
    }
    else if((freq > hi_pass) && (freq < hi_stop)) {
      v = 1.0 - (freq - hi_pass) / (hi_stop - hi_pass);
      PH[i] = v; 
    }
    else {
      // I have no idea... 
    }
  }

  FFT::shift(PH, PH);
}

// make the filter's impulse response.
// Frequencies may be negative.
// In fact, the filter's response will be symmetric around f = 0. 
SoDa::Filter::Filter(FTYPE typ, int num_taps, float sample_rate, float f1, float f2, 
		     float transition_width, float stopband_atten) :
  num_taps(num_taps) {

  float ny_lim = sample_rate / 2.0; 
  if((abs(f1) > ny_lim) || (abs(f2) > ny_lim)) {
  }

  // We *have* our standards, you know.
  if(stopband_atten < 50.0) stopband_atten = 50.0;
  
  // set all the relevant sizes
  h.resize(num_taps);

  H.resize(num_taps);
  // we haven't done an overlap and save op yet. 
  last_invec_size = 0;
  // we don't need an fft widget for the input stuff yet.
  dft_p = nullptr; 

  // make the prototype
  switch(typ) {
  case LO:
    makePrototype(H, sample_rate, 0.0, f1, transition_width);
    break;
  case HI:
    makePrototype(H, sample_rate, f1, 0.5 * sample_rate - transition_width, transition_width);    
    break;
  case BAND:
    makePrototype(H, sample_rate, f1, f2, transition_width);        
    break; 
  }

  // take the DFT
  SoDa::FFT fft(h.size());
  fft.ifft(h, H);

  // window it
  std::vector<float> cwin(num_taps);
  SoDa::ChebyshevWindow(cwin, num_taps, stopband_atten);

  // apply gain correction
  float gain_corr = 1.0 / ((float) num_taps); 
  for(int i = 0; i < num_taps; i++) {
    h[i] = h[i] * cwin[i] * gain_corr;
  }
}

void SoDa::Filter::makeImage(int image_size) {
  H.clear();

  H.resize(image_size);

  // create a new FFT object.
  if(dft_p != nullptr) {
    delete dft_p; 
  }

  std::cerr << "Making new FFT object sized " << image_size << "\n";
  dft_p = new SoDa::FFT(image_size); 

  std::vector<std::complex<float>> h_padded(image_size);

  int i, j; 
  // create the filter image.
  for(i = 0; i < image_size; i++) {
    h_padded[i] = std::complex<float>(0.0, 0.0);
  }

  for(i = 0; i < (num_taps + 1) / 2; i++, j++) {
    h_padded[i] = h[i]; 
  }
  for(i = 0; i < (num_taps / 2); i++) {
    h_padded[h_padded.size() - 1 - i] = h[h.size() - 1 - i];
  }

  temp.resize(image_size);
  
  // now transform it
  dft_p->fft(H, h_padded); 
}

void SoDa::Filter::setupNewOverlap(int invec_size) {
  debug_count = 0; 
  
  // how big should the transform be?
  // Make it a power of two or a power of two times 3 or 5.
  int target_min = invec_size + num_taps;
  int target_max = invec_size * 3 / 2; 

  // find a "good" size.
  // this is bone-brained stupid, but has the advantage that
  // I can explain it easily -- anything that's larger than
  // 16M entries isn't for us anyway.
  int good_size = 0;
  for(int i = 1; i < target_max; i = i * 2) {
    for(int j = 1; j < 8; j += 2) {
      if((i * j) > target_min) {
	int cand = i * j;
	if((cand < good_size) || (good_size == 0)) good_size = cand; 
      }
    }
  }

  if(good_size == 0) {
    std::stringstream ss;
    ss << "SoDa::Filter:\n\tCould not find a good overlap-save length for a " 
       <<  num_taps << " tap filter and " 
       << invec_size << "%1 entry vector\n";
    throw std::runtime_error(ss.str());
  }
  
  if((good_size - invec_size) > (invec_size / 3)) {
    std::stringstream ss;
    ss << "SoDa::Filter:\n\tGot a really bad overlap-save length of " 
       << good_size 
       << " for a " << num_taps  << " tap filter and " 
       << invec_size << " entry vector\n";
    throw std::runtime_error(ss.str());
  }

  

  // clear the saved buffer
  for(int i = 0; i < saved.size(); i++) {
    saved[i] = std::complex<float>(0.0, 0.0);
  }
  
  saved_dft.resize(good_size);
  
  save_length = good_size - invec_size; 
  
  last_invec_size = invec_size;

  saved.resize(good_size);
  
  makeImage(good_size);
}

  
