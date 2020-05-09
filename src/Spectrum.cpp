// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include "Spectrum.hpp"

using namespace uzume::dsp;

Spectrum::Spectrum(unsigned int fftSize) : spectralEnvelope(nullptr), aperiodicRatio(nullptr), fftSize(fftSize) {
    if (fftSize != 0) {
        spectralEnvelope = new double[fftSize];
        aperiodicRatio = new double[fftSize];
    }
}

Spectrum::~Spectrum() {
    delete[] spectralEnvelope;
    delete[] aperiodicRatio;
}
