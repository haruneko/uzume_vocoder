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

const std::function<AnalyzeAperiodicity *(unsigned int)> WaveformSpectrogram::DefaultAperiodicAnalysisFactory = [](unsigned int fftSize) {
    return new AnalyzeAperiodicityWithD4C(fftSize);
};

const std::function<AnalyzePeriodicity *(unsigned int)> WaveformSpectrogram::DefaultPeriodicAnalysisFactory = [](unsigned int fftSize) {
    return new AnalyzePeriodicityWithCheapTrick(fftSize);
};

const std::function<EstimateF0 *(double)> WaveformSpectrogram::DefaultF0EstimationFactory = [](double msFramePeriod) {
    return new EstimateF0WithDIO(msFramePeriod);
};

WaveformSpectrogram::WaveformSpectrogram(Waveform *waveform,
        const std::function<AnalyzeAperiodicity *(unsigned int)> &aperiodicAnalysisFactory,
        const std::function<AnalyzePeriodicity *(unsigned int)> &periodicAnalysisFactory,
        const std::function<EstimateF0 *(double)> &f0EstimationFactory)
        : analyzeAperiodicity(nullptr), analyzePeriodicity(nullptr), waveform(waveform), f0(nullptr), iw(nullptr) {

    analyzeAperiodicity = aperiodicAnalysisFactory(waveform->samplingFrequency);
    analyzePeriodicity = periodicAnalysisFactory(waveform->samplingFrequency);
    double msFramePeriod = 2.0;

    auto f0Estimation = f0EstimationFactory(msFramePeriod);
    f0 = new Contour(waveform->msLength(), msFramePeriod);
    (*f0Estimation)(f0, waveform);
    delete f0Estimation;

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
