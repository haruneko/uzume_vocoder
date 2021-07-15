// Copyright 2021 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include "ArraySpectrogramAggregator.hpp"

using namespace uzume::vocoder;

ArraySpectrogramAggregator::ArraySpectrogramAggregator(std::vector<const Spectrogram *> &spectrograms)
    : spectrograms(spectrograms) {

}

ArraySpectrogramAggregator *ArraySpectrogramAggregator::from(std::vector<const Spectrogram *> &spectrograms) {
    if(spectrograms.empty()) {
        return nullptr;
    }
    for(auto s: spectrograms) {
        if(s->fftSize() != spectrograms[0]->fftSize()) {
            return nullptr;
        }
    }
    return new ArraySpectrogramAggregator(spectrograms);
}

std::pair<int, double> ArraySpectrogramAggregator::indexAndMsAt(double ms) const {
    double currentMs = 0.0;
    double previousMs;
    for(int i = 0; i < spectrograms.size(); i++) {
        previousMs = currentMs;
        currentMs += spectrograms[i]->msLength();
        if(ms <= currentMs) {
            return std::pair<int, double>(i, ms - previousMs);
        }
    }
    return std::pair<int, double>(spectrograms.size() - 1, currentMs);
}

double ArraySpectrogramAggregator::f0At(double ms) const {
    if(spectrograms.empty()) {
        return 0.0;
    }

    auto p = indexAndMsAt(ms < 0.0 ? 0.0 : ms);
    return spectrograms[p.first]->f0At(p.second);
}

bool ArraySpectrogramAggregator::pickUpSpectrumAt(Spectrum *destination, double ms) const {
    if(spectrograms.empty()) {
        return false;
    }
    auto p = indexAndMsAt(ms);
    return spectrograms[p.first]->pickUpSpectrumAt(destination, p.second);
}

double ArraySpectrogramAggregator::msLength() const {
    double r = 0.0;
    for(auto s: spectrograms) {
        r += s->msLength();
    }
    return r;
}

unsigned int ArraySpectrogramAggregator::fftSize() const {
    if(spectrograms.empty()) {
        return 0;
    }
    return spectrograms[0]->fftSize();
}
