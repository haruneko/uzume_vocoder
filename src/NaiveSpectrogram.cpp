// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include <cmath>
#include "NaiveSpectrogram.hpp"
#include "util.hpp"

using namespace uzume::dsp;

namespace {
    inline void initializeSpecgram(double ***target, unsigned int length, unsigned int fftSize) {
        *target = new double *[length];
        for (int i = 0; i < length; i++) {
            (*target)[i] = new double[fftSize / 2 + 1];
        }
    }

    inline void deleteSpecgram(unsigned int length, double ***target) {
        for(int i = 0; i < length; i++) {
            delete[] (*target)[i];
        }
        delete[] *target;
    }

    template<class T>
    inline T myMin(T x, T y) {
        return (x < y) ? x : y;
    }

    template<class T>
    inline T myMax(T x, T y) {
        return (x > y) ? x : y;
    }

    inline double safeAperiodicity(double x) {
        return myMax<double>(0.001, myMin<double>(0.999999999999, x));
    }
} // namespace

NaiveSpectrogram::NaiveSpectrogram(unsigned int length, unsigned int fftSize, double msFramePeriod)
        : periodicSpecgram(nullptr), aperiodicSpecgram(nullptr), f0Contour(nullptr), timeAxis(nullptr),
          length(length), _fftSize(fftSize), msFramePeriod(msFramePeriod) {
    if (length != 0 && fftSize != 0) {
        initializeSpecgram(&periodicSpecgram, length, fftSize);
        initializeSpecgram(&aperiodicSpecgram, length, fftSize);
        f0Contour = new double[fftSize];
        timeAxis = new double[fftSize];
    }
}

NaiveSpectrogram::~NaiveSpectrogram() {
    deleteSpecgram(length, &periodicSpecgram);
    deleteSpecgram(length, &aperiodicSpecgram);
    delete[] f0Contour;
    delete[] timeAxis;
}

bool NaiveSpectrogram::pickUpSpectrumAt(Spectrum *destination, double ms) const {
    if (destination->fftSize != _fftSize) {
        return false;
    }
    int currentFrameIndexFloor = myMin<int>((int) (floor(ms / msFramePeriod)), length - 1);
    int currentFrameIndexCeil = myMin<int>((int) (ceil(ms / msFramePeriod)), length - 1);
    double interpolation = ms / msFramePeriod - currentFrameIndexFloor;
    for (unsigned int i = 0; i <= _fftSize / 2; i++) {
        destination->periodicSpectrum[i] =
                fabs(periodicSpecgram[currentFrameIndexFloor][i] * (1.0 - interpolation)) +
                fabs(periodicSpecgram[currentFrameIndexCeil][i] * interpolation);
    }
    for (unsigned int i = 0; i <= _fftSize / 2; i++) {
        destination->aperiodicSpectrum[i] =
                pow(safeAperiodicity(aperiodicSpecgram[currentFrameIndexFloor][i]) * (1.0 - interpolation) +
                    safeAperiodicity(aperiodicSpecgram[currentFrameIndexCeil][i]) * interpolation, 2.0);
    }
    return true;
}

double NaiveSpectrogram::f0At(double ms) const {
    auto f0IndexFloor = (int)(floor(ms / msFramePeriod));
    auto f0IndexCeil = (int)(ceil(ms / msFramePeriod));
    double interpolation = (ms - f0IndexFloor * msFramePeriod) / msFramePeriod;
    return f0Contour[f0IndexFloor] * (1.0 - interpolation) + f0Contour[f0IndexCeil] * interpolation;
}

double NaiveSpectrogram::msLength() const {
    return length * msFramePeriod;
}

unsigned int NaiveSpectrogram::fftSize() const {
    return _fftSize;
}
