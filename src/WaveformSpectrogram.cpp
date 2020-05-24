// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include <algorithm>
#include "AnalyzeAperiodicityWithD4C.hpp"
#include "AnalyzePeriodicity.hpp"
#include "EstimateF0WithDIO.hpp"
#include "util.hpp"
#include "WaveformSpectrogram.hpp"
#include "AnalyzePeriodicityWithCheapTrick.hpp"

using namespace uzume::dsp;

WaveformSpectrogram::WaveformSpectrogram(Waveform *waveform)
        : analyzeAperiodicity(nullptr), analyzePeriodicity(nullptr), waveform(waveform), f0(nullptr), iw(nullptr) {
    // FIXME: needs DI to new those analysis classes.
    analyzeAperiodicity = new AnalyzeAperiodicityWithD4C(waveform->samplingFrequency);
    analyzePeriodicity = new AnalyzePeriodicityWithCheapTrick(waveform->samplingFrequency);
    double msFramePeriod = 2.0;

    EstimateF0WithDIO dio(msFramePeriod);
    f0 = new Contour(waveform->msLength(), msFramePeriod);
    dio(f0, waveform);

    iw = new InstantWaveform(waveform->samplingFrequency);
}

WaveformSpectrogram::~WaveformSpectrogram() noexcept {
    delete analyzeAperiodicity;
    delete analyzePeriodicity;
    delete waveform;
    delete f0;
    delete iw;
}

unsigned int WaveformSpectrogram::fftSize() const {
    return uzume::dsp::fftSize(waveform->samplingFrequency);
}

double WaveformSpectrogram::msLength() const {
    return f0->msLength();
}

double WaveformSpectrogram::f0At(double ms) const {
    return f0->at(ms);
}

bool WaveformSpectrogram::pickUpSpectrumAt(Spectrum *destination, double ms) const {
    int origin = (int) (ms / 1000.0 * waveform->samplingFrequency);
    for (int i = -(int) fftSize(); i < (int) fftSize(); i++) {
        int waveIndex = std::max(0, std::min((int) waveform->length - 1, origin + i));
        iw->data[i] = waveform->data[waveIndex];
    }
    iw->f0 = f0At(ms);
    return (*analyzePeriodicity)(destination, iw) && (*analyzeAperiodicity)(destination, iw);
}
