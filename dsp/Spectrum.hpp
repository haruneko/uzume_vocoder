// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef __DSP_SPECTRUM_HPP
#define __DSP_SPECTRUM_HPP

/**
 * Spectrum represents a pair of spectra, periodic one and aperiodic one.
 */
class Spectrum {
public:
    Spectrum(unsigned int fftSize);

    ~Spectrum();

    double *spectralEnvelope;
    double *aperiodicRatio;
    unsigned int fftSize;
};

#endif
