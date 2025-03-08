/* this file exists primarily to create a mainpage.  */

/**
 */


#include "Filter.hxx"
#include "ReSampler.hxx"
#include "FFT.hxx"

///
/// @file Signals.hxx
///
/// @brief A Library of basic signal processing functions
///
/// @author M. H. Reilly (kb1vc)
/// @date Jan 2021,2025
///

/**
 * @mainpage SoDa Signals: A Library of basic signal processing functions
 * 
 * @section Overview
 * 
 * SoDa Signals provides simple classes and functions that perform FFT/IFFT, 
 * and apply filters to complex (I and Q) signals.  
 * 
 * All (or at least *almost all*) of the functions and classes here are 
 * aimed at manipulating quadrature signals with an in-phase (I) and a
 * quadrature (Q) component.
 *
 * If quadrature signal representation is unfamiliar, 
 * it would be to your advantage to learn about it.
 * The wikipedia article on 
 * "In-phase and quadrature components" is particularly inadequate, in fact
 * it is worth ignoring.
 *
 * Lyons offers a good presentation in "Understanding
 * Digital Signal Processing."
 * In fact, it is so good, that one wonders why
 * he waited until chapter eight to introduce it. 
 *
 * In any case, almost all the interesting modulation and demodulation schemes, 
 * filtering, and frequency translation, work better if the signal is represented
 * as a complex number. 
 *
 * This library is used in the SoDaRadio SDR program. So it spends
 * a whole lot of effort on real-time performance. SoDaRadio runs on
 * some pretty primitive hardware: keeping compute costs down is
 * important. 
 * 
 * 
 * @section FFT The SoDa::FFT Class
 * 
 * SoDa::FFT is a sanitized interface to the FFTW (Fastest Fourier Transform in the West)
 * library.  The library offers pretty good performance, though perhaps not quite as
 * good as the Intel MKL implementations. When this library was written, the Intel MKL
 * library did not have a good anwer for non-power-of-two DFTs. Unfortunately, the OSFilter
 * and ReSampler classes really need odd size DFTs to make real time operations on
 * SDR sample streams work with standard audio rates. (Implementing a 1.25 MHz to 48 kHz
 * resampler without a prime-radix DFT gets costly quick.)
 *
 * The SoDa::FFT class provides transforms for both complex<double> and complex<float>
 * signals. It handles much of the boilerplate and fiddly bits around things like
 * buffer allocation and such. It even handles training and tuning. 
 *
 * @section Filter The SoDa::Filter Class
 *
 * Almost any DSP chain will need a filter somewhere along the line. Textbooks
 * are replete with discussions of IIR and FIR filters, Z transforms, Parks-McClellan, 
 * and all that stuff.  
 *
 * Anybody attempting to learn digital signal processing should
 * implement a few filters on their own.  However, the SoDa::Filter
 * class can provide a reference or provide a check on a student's
 * experiment. 
 *
 * SoDa::Filter and SoDa::OSFilter (and SoDa::Resampler) all perform
 * their operations in the frequency domain.  It turns out that for any
 * filter larger than about 8 taps, doing the processing with a DFT
 * is far more efficient. By a lot. 
 *
 * SoDa::Filter processes just one input block at a time, so it isn't
 * directly useful for streaming data. If you try to use it that way,
 * there will be unpleasant artifacts at the start of each output buffer.
 * Look it up -- the problem has to do with circular convolution that
 * came along with our DFT based filter scheme. 
 *
 * SoDa::Filter constructs a specified filter using the Window method. This is
 * much more flexible than Parks-McClellan or Remez schemes. They were good at
 * minimizing the number of taps required for a filter, but give that we're doing
 * DFT based filters on reasonably sized buffers (> 1K samples is good) the
 * number of taps required for a filter is nearly irrelevant. And neither of
 * those methods is very good at strange filter shapes. 
 *
 * @section OSFilter The SoDa::OSFilter Class - so much more useful than SoDa::Filter
 *
 * SoDa::OSFilter operates on continuous streams --
 * using the overlap-and-save method. (Look at the wikipedia page and marvel at how
 * much it overthinks the whole thing.)
 *
 * Remember that "circular convolution" thing that was mentioned in @ref Filter? You
 * should already be familiar with the concept of an FIR filter.  If not, do some
 * reading and come back to this later.
 *
 * An FIR filter applies its M filter taps to a vector constructed of the last M samples
 * in the input stream.  But that buffer is empty for the first M-1 samples.  If you're
 * implementing your filter in the time domain on a single block (and why should you?) you'd ignore the
 * first M-1 samples and say the heck with it.  But doing that for a continuous stream
 * causes nasty bursts of noise at the start of each buffer. So, we can avoid all that
 * by "priming the pump." We load the last M-1 samples into the filter's input buffer
 * and *then* throw away the first M-1 output samples.  But the rest of the stuff will
 * be just fine.
 *
 * That's all a consequence of "circular convolution."  And we can't escape it simply
 * by doing everything with DFTs. The OSFilter class processes streams of blocks, saving
 * the tail end of each block to prepend at the start of the next block.
 *
 * 
 * @section ReSampling The SoDa::ReSampler Class
 *
 * In a software defined radio, samples come in at a rate that is
 * determined by the radio hardware at the other end, bandwidth over the wire,
 * channel requirements, and all that jazz.  Sometimes that rate isn't
 * what we need. Few radios support sample rates of 48 kS/s -- a
 * fairly common sample rate for audio widgets. We may even want to
 * do some processing before writing a signal stream to a disk or
 * sending it out the audio port.
 * 
 * That's where the ReSampler class comes in.
 *
 * Like SoDa::Filter, SoDa::ReSampler operates on a continuous signal
 * stream. Otherwise it would be pretty useless. 
 *
 * @section SoDa
 * 
 * SoDa is a namespace around a set of classes, libraries, (and one
 * application) that have been developed as a hobby.  They are not
 * meant for production (though the license is permissive here) as
 * they are a running chronicle of my own -- sometimes slow -- process
 * of learning the knooks and krannies of DSP.
 * 
 * SoDaSignals is used in SoDaRadio.
 * 
 * SoDaRadio is the first real DSP application I wrote.  It is a
 * software-defined radio application written for the USRP line of
 * SDRs from Ettus Research. I've used it as an "all-mode" amateur
 * radio on bands ranging from 10 MHz to 10 GHz.  Much of the code in
 * SoDaSignals will eventually replace components in SoDaRadio.
 * (Especially the really grody and barely competent implementation of
 * the filters.)
 *
 * All SoDa components are FOSS.  See the associated licenses.
 *
 * @section License
 *
 *  @copyright  2025 Matt Reilly - kb1vc  All rights reserved.
 *
 *  @par BSD 2-Clause License
 *
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice, this
 *  list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
