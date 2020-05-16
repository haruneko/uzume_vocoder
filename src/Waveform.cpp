// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include "Waveform.hpp"

using namespace uzume::dsp;

Waveform::Waveform(unsigned int length, unsigned int samplingFrequency)
        : data(nullptr), length(length), samplingFrequency(samplingFrequency){
    if(length > 0) {
        data = new double[length];
    }
}

Waveform::~Waveform() {
    delete[] data;
}
