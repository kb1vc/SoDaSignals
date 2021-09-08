/**
 * @file SoDaSignals.i
 * 
 * @brief swig wrapper for SoDaSignals classes
 */

/* start with FFT class */



%module SoDaSignals
%{
#define SWIG_FILE_WITH_INIT
  /* #include "FFT.hxx" */
#include "Signals.hxx"
%}

%include <std_vector.i>
%include <std_complex.i>

namespace std {
  %template(vectorcf) vector<complex<float>>;
  %template(vectorcd) vector<complex<double>>;
  %template(vectorf) vector<float>;
  %template(vectord) vector<double>;
 };

namespace SoDa {
  %template(FilterCF) Filter<float>;
  %template(FilterCD) Filter<double>;
};


%include "FFT.hxx" 
%include "Signals.hxx"
%include "Filter.hxx"

