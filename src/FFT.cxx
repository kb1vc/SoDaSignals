#include "FFT.hxx"
SoDa::FFT::FFT(size_t N, unsigned int flags) {
  dim = N;

  fftw_set_timelimit(1.0);
  fftwf_set_timelimit(1.0);

  f_dummy_in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * dim);
  f_dummy_out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * dim);  
  d_dummy_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim);
  d_dummy_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim);  
  
  fplan_f = fftwf_plan_dft_1d(dim, f_dummy_in, f_dummy_out, FFTW_FORWARD, flags);

  fplan_i = fftwf_plan_dft_1d(dim, f_dummy_in, f_dummy_out, FFTW_BACKWARD, flags);

  dplan_f = fftw_plan_dft_1d(dim, d_dummy_in, d_dummy_out, FFTW_FORWARD, flags);

  dplan_i = fftw_plan_dft_1d(dim, d_dummy_in, d_dummy_out, FFTW_BACKWARD, flags);
  
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
  fftwf_execute_dft(fplan_f, (fftwf_complex *) in.data(), (fftwf_complex *) out.data());
  return true;
}

bool SoDa::FFT::ifft(std::vector<std::complex<float>> & out, std::vector<std::complex<float>> & in)
{
  fftwf_execute_dft(fplan_i, (fftwf_complex *) in.data(), (fftwf_complex *) out.data());
  return true;
}



bool SoDa::FFT::fft(std::vector<std::complex<double>> & out, std::vector<std::complex<double>> & in)
{
  fftw_execute_dft(dplan_f, (fftw_complex *) in.data(), (fftw_complex *) out.data());
  return true;
}

bool SoDa::FFT::ifft(std::vector<std::complex<double>> & out, std::vector<std::complex<double>> & in)
{
  fftw_execute_dft(dplan_i, (fftw_complex *) in.data(), (fftw_complex *) out.data());
  return true;
}


