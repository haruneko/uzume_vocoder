// Copyright 2021 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include "StretchedPartialSpectrogram.hpp"

using namespace uzume::vocoder;

StretchedPartialSpectrogram::StretchedPartialSpectrogram(const Spectrogram *spectrogram, const TimeAxisMap *timeAxisMap)
    : spectrogram(spectrogram), timeAxisMap(timeAxisMap) {
}

bool StretchedPartialSpectrogram::pickUpSpectrumAt(Spectrum *destination, double ms) const {
    return spectrogram->pickUpSpectrumAt(destination, timeAxisMap->at(ms));
}

double StretchedPartialSpectrogram::f0At(double ms) const {
    return spectrogram->f0At(timeAxisMap->at(ms));
}

double StretchedPartialSpectrogram::msLength() const {
    return timeAxisMap->msLength();
}

unsigned int StretchedPartialSpectrogram::fftSize() const {
    return spectrogram->fftSize();
}
