/**
 * @file SoDaSignals.i
 * 
 * @brief swig wrapper for SoDaSignals classes
 */

/* start with FFT class */



%module SoDaSignals
%{
#include "FFT.hxx"
%}

%include <std_vector.i>
%include <std_complex.i>
namespace std {
  %template(vectorcf) vector<complex<float>>;
  %template(vectorcd) vector<complex<double>>;
  %template(vectorf) vector<float>;
  %template(vectord) vector<double>;

};

%include "FFT.hxx"
