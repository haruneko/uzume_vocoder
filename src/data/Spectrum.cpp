// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include "Spectrum.hpp"

using namespace uzume::vocoder;

Spectrum::Spectrum(unsigned int fftSize) : periodicSpectrum(nullptr), aperiodicSpectrum(nullptr), fftSize(fftSize) {
    if (fftSize != 0) {
        periodicSpectrum = new double[fftSize];
        aperiodicSpectrum = new double[fftSize];
    }
}

Spectrum::~Spectrum() {
    delete[] periodicSpectrum;
    delete[] aperiodicSpectrum;
}
