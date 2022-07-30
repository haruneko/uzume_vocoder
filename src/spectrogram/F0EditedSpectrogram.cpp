// Copyright 2022 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include <cmath>
#include "F0EditedSpectrogram.hpp"

using namespace uzume::vocoder;

const std::function<double(double, double)> linear = [](double f0, double ratio) {
    return f0 * ratio;
};

const std::function<double(double, double)> logBy2 = [](double f0, double ratio) {
    return pow(2, ratio) * f0;
};

F0EditedSpectrogram::F0EditedSpectrogram(const Spectrogram *spectrogram, const ControlChange &cc, const SynthType synthType)
: spectrogram(spectrogram), cc(cc) {
    switch(synthType) {
    case SynthType::Linear:
        this->f = linear;
        break;
    case SynthType::Log:
        this->f = logBy2;
        break;
    }
}

bool F0EditedSpectrogram::pickUpSpectrumAt(Spectrum *destination, double ms) const {
    return spectrogram->pickUpSpectrumAt(destination, ms);
}

double F0EditedSpectrogram::f0At(double ms) const {
    double position = ms / msLength();
    return this->f(cc.at(position), spectrogram->f0At(ms));
}

double F0EditedSpectrogram::msLength() const {
    return spectrogram->msLength();
}

unsigned int F0EditedSpectrogram::fftSize() const {
    return spectrogram->fftSize();
}
