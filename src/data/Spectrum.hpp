// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef __DSP_SPECTRUM_HPP
#define __DSP_SPECTRUM_HPP

namespace uzume { namespace dsp {
/**
 * Spectrum represents a pair of spectra, periodic one and aperiodic one.
 * This class HAS a responsibility to free its allocated memories.
 */
class Spectrum final {
public:
    Spectrum(unsigned int fftSize);

    ~Spectrum();

    double *periodicSpectrum;
    double *aperiodicSpectrum;
    unsigned int fftSize;
};

} }

#endif
