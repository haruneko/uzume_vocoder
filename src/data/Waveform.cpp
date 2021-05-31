// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include "../world/audioio.h"
#include "Waveform.hpp"

using namespace uzume::vocoder;

Waveform::Waveform(unsigned int length, unsigned int samplingFrequency)
        : data(nullptr), length(length), samplingFrequency(samplingFrequency){
    if(length > 0) {
        data = new double[length];
        this->clear();
    }
}

Waveform::~Waveform() {
    delete[] data;
}

double Waveform::msLength() const {
    return (double)length / samplingFrequency * 1000.0;
}

Waveform *Waveform::read(const char *filepath) {
    int waveLength = GetAudioLength(filepath);
    auto *res = new Waveform(waveLength, 44100);
    int fs, nbits;
    wavread(filepath, &fs, &nbits, res->data);
    res->samplingFrequency = (unsigned int)fs;
    return res;
}

bool Waveform::save(const char *filepath, int bits) const {
    wavwrite(data, (int)length, (int)samplingFrequency, bits, filepath);
    return true;
}

void Waveform::clear() {
    if(data) {
       for(unsigned int i = 0; i < length; i++) {
           data[i] = 0.0;
       }
    }
}
