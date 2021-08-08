// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include "../world/audioio.h"
#include "Waveform.hpp"
#include <cmath>

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

Waveform *Waveform::read(const char *filepath) {
    int waveLength = GetAudioLength(filepath);
    if(waveLength <= 0) {
        return nullptr;
    }
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

double Waveform::maxAbsoluteValueBetween(double msBegin, double msEnd) const {
    return maxAbsoluteValueBetween(indexAt(msBegin), indexAt(msEnd));
}

double Waveform::maxAbsoluteValueBetween(int indexBegin, int indexEnd) const {
    double r = 0.0;
    int indexStart = indexBegin < 0 ? 0 : indexBegin;
    int indexStop = indexEnd < length ? indexEnd : (int)length;
    for(int i = indexStart; i < indexStop; i++) {
        double v = fabs(data[i]);
        r = r < v ? v : r;
    }
    return r;
}

double Waveform::rootMeanSquareBetween(double msBegin, double msEnd) const {
    return rootMeanSquareBetween(indexAt(msBegin), indexAt(msEnd));
}

double Waveform::rootMeanSquareBetween(int indexBegin, int indexEnd) const {
    double r = 0.0;
    int indexStart = indexBegin < 0 ? 0 : indexBegin;
    int indexStop = indexEnd < length ? indexEnd : (int)length;
    for(int i = indexStart; i < indexStop; i++) {
        r += data[i] * data[i];
    }
    if(indexStop - indexStart != 0) {
        r /= (double)(indexStop - indexStart);
    }
    return sqrt(r);
}
