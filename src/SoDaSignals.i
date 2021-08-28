/**
 * @file SoDaSignals.i
 * 
 * @brief swig wrapper for SoDaSignals classes
 */

/* start with FFT class */

/* %include <std_vector.i> */

%module SoDaSignals
%{
#include "FFT.hxx"
%}

%include "FFT.hxx"
