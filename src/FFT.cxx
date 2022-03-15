#include "FFT.hxx"
#include <sstream>
#include <iostream>
#include <exception>
#include <stdexcept>

SoDa::FFT::FFT(size_t N, unsigned int flags,
	       int istride, 
	       int ostride) {
  dim = N;

  w_float.resize(N);
  w_double.resize(N);
  s_cfloat.resize(N);
  s_cdouble.resize(N);
  
  fftw_set_timelimit(1.0);
  fftwf_set_timelimit(1.0);

  f_dummy_in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * dim * istride);
  f_dummy_out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * dim * ostride);  
  d_dummy_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim * istride);
  d_dummy_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim * ostride);  

  if((istride == 1) && (ostride == 1)) {
    fplan_f = fftwf_plan_dft_1d(dim, f_dummy_in, f_dummy_out, FFTW_FORWARD, flags);

    fplan_i = fftwf_plan_dft_1d(dim, f_dummy_in, f_dummy_out, FFTW_BACKWARD, flags);

    dplan_f = fftw_plan_dft_1d(dim, d_dummy_in, d_dummy_out, FFTW_FORWARD, flags);

    dplan_i = fftw_plan_dft_1d(dim, d_dummy_in, d_dummy_out, FFTW_BACKWARD, flags);
  }
  else {
    int np = dim; 
    int * npp = & np; 
    // we are planning several parallel DFTs.
    fplan_f = fftwf_plan_many_dft(1, npp, istride, 
				  f_dummy_in, npp,
				  istride, 1, 
				  f_dummy_out, npp,
				  ostride, 1,
				  FFTW_FORWARD, 				  
				  flags);
    fplan_i = fftwf_plan_many_dft(1, npp, istride, 
				  f_dummy_in, npp,
				  istride, 1, 
				  f_dummy_out, npp,
				  ostride, 1,
				  FFTW_BACKWARD,
				  flags);

    dplan_f = fftw_plan_many_dft(1, npp, istride, 
				  d_dummy_in, npp,
				  istride, 1, 
				  d_dummy_out, npp,
				  ostride, 1,
				  FFTW_FORWARD, 				  
				  flags);
    dplan_i = fftw_plan_many_dft(1, npp, istride, 
				  d_dummy_in, npp,
				  istride, 1, 
				  d_dummy_out, npp,
				  ostride, 1,
				  FFTW_BACKWARD,
				  flags);
  }
  
}

SoDa::FFT::~FFT() {
  fftwf_destroy_plan(fplan_f);
  fftwf_destroy_plan(fplan_i);  
  fftw_destroy_plan(dplan_f);
  fftw_destroy_plan(dplan_i);  
  
  fftwf_free(f_dummy_in);
  fftwf_free(f_dummy_out);  
  fftw_free(d_dummy_in);
  fftw_free(d_dummy_out);  
}

bool SoDa::FFT::fft(std::vector<std::complex<float>> & out, std::vector<std::complex<float>> & in)
{
  checkInOut(out.size(), in.size());
  fftwf_execute_dft(fplan_f, (fftwf_complex *) in.data(), (fftwf_complex *) out.data());
  return true;
}

bool SoDa::FFT::ifft(std::vector<std::complex<float>> & out, std::vector<std::complex<float>> & in)
{
  checkInOut(out.size(), in.size());  
  fftwf_execute_dft(fplan_i, (fftwf_complex *) in.data(), (fftwf_complex *) out.data());
  return true;
}



bool SoDa::FFT::fft(std::vector<std::complex<double>> & out, std::vector<std::complex<double>> & in)
{
  checkInOut(out.size(), in.size());
  fftw_execute_dft(dplan_f, (fftw_complex *) in.data(), (fftw_complex *) out.data());
  return true;
}

bool SoDa::FFT::ifft(std::vector<std::complex<double>> & out, std::vector<std::complex<double>> & in)
{
  checkInOut(out.size(), in.size());  
  fftw_execute_dft(dplan_i, (fftw_complex *) in.data(), (fftw_complex *) out.data());
  return true;
}

bool SoDa::FFT::spectrogram(std::vector<std::complex<double>> & out, std::vector<std::complex<double>> & in)
{
  checkInOut(out.size(), in.size());  
  // apply the window
  for(int i = 0; i < in.size(); i++) {
    s_cdouble[i] = in[i] * w_double[i];
  }
  fftw_execute_dft(dplan_f, (fftw_complex *) s_cdouble.data(), (fftw_complex *) out.data());
  
  return true;
}

bool SoDa::FFT::spectrogram(std::vector<std::complex<float>> & out, std::vector<std::complex<float>> & in)
{
  checkInOut(out.size(), in.size());  
  // apply the window
  for(int i = 0; i < in.size(); i++) {
    s_cfloat[i] = in[i] * w_float[i];
  }
  fftwf_execute_dft(fplan_f, (fftwf_complex *) s_cfloat.data(), (fftwf_complex *) out.data());
  
  return true;
}


void SoDa::FFT::initBlackmanHarris(int size)
{
  unsigned int i;
  double a0 = 0.35875;
  double a1 = 0.48829;
  double a2 = 0.14128;
  double a3 = 0.01168;

  w_float.resize(size);
  w_double.resize(size);
  float anginc = 2.0 * M_PI / ((float) size - 1);
  float ang = 0.0; 
  for(i = 0; i < size; i++) {
    ang = anginc * ((float) i); 
    w_float[i] = a0 - a1 * cos(ang) + a2 * cos(2.0 * ang) -a3 * cos(3.0 * ang);
    w_double[i] = a0 - a1 * cos(ang) + a2 * cos(2.0 * ang) -a3 * cos(3.0 * ang);     
  }
}



void SoDa::FFT::checkInOut(size_t outsize, size_t insize) {
  std::stringstream ss; 
  if((outsize == dim) && (insize == dim)) return;
  
  if(outsize != dim) {
    ss << "Output vector size " << outsize << " "; 
  }
  if(insize != dim) {
    ss << "Input vector size " << insize << " ";
  }
  ss << "is not equal to the planned size of " << dim << " elements\n";
  throw std::runtime_error(ss.str());
}
