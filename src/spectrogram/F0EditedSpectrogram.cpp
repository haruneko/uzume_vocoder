// Copyright 2022 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include "F0EditedSpectrogram.hpp"

using namespace uzume::vocoder;

F0EditedSpectrogram::F0EditedSpectrogram(const Spectrogram *spectrogram, const ControlChange &cc)
    : spectrogram(spectrogram), cc(cc) {
}

bool F0EditedSpectrogram::pickUpSpectrumAt(Spectrum *destination, double ms) const {
    return spectrogram->pickUpSpectrumAt(destination, ms);
}

double F0EditedSpectrogram::f0At(double ms) const {
    double position = ms / msLength();
    return cc.at(position) * spectrogram->f0At(ms);
}

double F0EditedSpectrogram::msLength() const {
    return spectrogram->msLength();
}

unsigned int F0EditedSpectrogram::fftSize() const {
    return spectrogram->fftSize();
}
