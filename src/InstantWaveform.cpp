// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include "InstantWaveform.hpp"

#include "util.hpp"

using namespace uzume::dsp;

InstantWaveform::InstantWaveform(unsigned int samplingFrequency)
        : data(nullptr), length(0), samplingFrequency(samplingFrequency) {
    length = uzume::dsp::fftSize(samplingFrequency);
    data = (new double[length * 2]) + length;
    f0 = 0.0;
}

InstantWaveform::~InstantWaveform() {
    delete[] (data - length);
}
