// Copyright 2021 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.

#include "data/Waveform.hpp"
#include "spectrogram//StretchedPartialSpectrogram.hpp"
#include "spectrogram/WaveformSpectrogram.hpp"

using namespace uzume::vocoder;

int main(void) {
    auto *in = Waveform::read("hogeeee");
    auto *inSpec = new WaveformSpectrogram(in);

    class Tam : public TimeAxisMap {
        double at(double ms) const override final {
            return ms * 2.0;
        };
        double msLength() const override final {
            return 1000.0;
        }
    };
    auto *tam = new Tam();
    auto *outSpec = new StretchedPartialSpectrogram(inSpec, tam);
    auto *out = new Waveform(in->length * 2, in->samplingFrequency);
    return 0;
}