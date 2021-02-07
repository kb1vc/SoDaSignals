/* this file exists solely to create a mainpage.  */

///
/// @file Signals.hxx
///
/// @brief A Library of basic signal processing functions
///
/// @author M. H. Reilly (kb1vc)
/// @date Jan 2021
///

/**
 * @mainpage SoDa Signals: A Library of basic signal processing functions
 * 
 * \section Overview
 * 
 * SoDa Signals provides simple classes and functions that perform FFT/IFFT, 
 * and apply filters to complex (I and Q) signals.  
 * 
 * All (or at least *almost all*) of the functions and classes here are 
 * aimed at manipulating quadrature signals with an in-phase (I) and a
 * quadrature (Q) component.  If quadrature signal representation is unfamiliar, 
 * it would be to your advantage to learn about it. The wikipedia article on 
 * "In-phase and quadrature components" is particularly inadequate, in fact
 * it is worth ignoring.  Lyons offers a good presentation in "Understanding
 * "Digital Signal Processing."  In fact, it is so good, that one wonders why
 * he waited until chapter eight to introduce it. 
 * 
 * In any case, almost all the interesting modulation and demodulation schemes, 
 * filtering, and frequency translation, work better if the signal is represented
 * as a complex number. 
 * 
 * \section FFT The SoDa::FFT Class
 * 
 * SoDa::FFT is a sanitized interface to the FFTW (Fastest Fourier Transform in the West)
 * library.  The library offers pretty good performance, though perhaps not quite as
 * good as the Intel MKL implementations. But Intel doesn't do a very good job of 
 * making MKL easy to install (as of Jan 2021) or even find. 
 * 
 * So SoDa Signals uses FFTW3. 
 * 
 * The SoDa::FFT class provides transforms for both complex<double> and complex<float>
 * signals. 
 *
 * \section Filter The SoDa::Filter Class
 *
 * Almost any DSP chain will need a filter somewhere along the line. Textbooks
 * are replete with discussions of IIR and FIR filters, Z transforms, Parks-McClellan, 
 * and all that stuff.  
 *
 * Anybody attempting to learn digital signal processing should implement 
 * a few filters on their own.  However, the SoDa::Filter class can provide
 * a reference or provide a check on a student's experiment. 
 * 
 * \section Windows 
 * 
 * At this writing, there is but one window function in the library:
 * The Dolph-Chebyshev window.  This was seen as a good choice for a
 * filter-synthesis window (See Lyons "Understading Digital Signal
 * Processing" for a discussion of the "Window Design Method."  For a
 * definitive catalog of windows, see the paper by Fred Harris "On the
 * Use of Windows for Harmonic Analysis With the Discrete Fourier
 * Transform.")
 * 
 * The Dolph-Chebyshev window, is somewhat difficult to implement from
 * the common literature. Both of the above sources were incomplete, especially
 * regarding the definition of the Chebyshev function that is crucial
 * of the window's construction.  See the comments in the sources for 
 * ChebyshevWindow.cxx to get more details. 
 *
 * \section SoDa
 * 
 * SoDa is a set of classes, libraries, (and one application) that have
 * been developed as a hobby.  They are not meant for production (though
 * the license is permissive here) as they are a running chronicle of 
 * my own -- sometimes slow -- process of learning the knooks and krannies
 * of DSP. 
 * 
 * SoDaSignals is this collection of basic blocks. 
 * 
 * SoDaRadio is the first real DSP application I wrote.  It is a software-defined
 * radio application written for the USRP line of SDRs from Ettus Research. I've
 * used it as an "all-mode" amateur radio on bands ranging from 10 MHz to 10 GHz. 
 * Much of the code in SoDaSignals will eventually replace components in SoDaRadio. 
 * (Especially the really grody and barely competent implementation of the filters.)
 *
 * SoDa::Format is more generally useful: it is a package for constructing 
 * formatted strings that show floating point numbers the way God intended: in 
 * engineering notation where the power-of-ten exponent is a multiple of 3. 
 * 
 * All SoDa components are FOSS.  See the associated licenses.
 * 
 * 
 */
