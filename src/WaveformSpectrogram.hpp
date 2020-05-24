// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_DSP_WAVEFORMSPECTROGRAM_HPP
#define UZUME_DSP_WAVEFORMSPECTROGRAM_HPP
#include "AnalyzeAperiodicity.hpp"
#include "AnalyzePeriodicity.hpp"
#include "Contour.hpp"
#include "Spectrogram.hpp"
#include "Waveform.hpp"

namespace uzume { namespace dsp {

class WaveformSpectrogram final : public Spectrogram {
public:
    explicit WaveformSpectrogram(Waveform *waveform);
    ~WaveformSpectrogram() final;

    bool pickUpSpectrumAt(Spectrum *destination, double ms) const override ;

    double f0At(double ms) const override ;

    double msLength() const override ;

    unsigned int fftSize() const override ;

private:
    AnalyzeAperiodicity *analyzeAperiodicity;
    AnalyzePeriodicity *analyzePeriodicity;
    Waveform *waveform;
    Contour *f0;

    InstantWaveform *iw;
};

} }

#endif //UZUME_DSP_WAVEFORMSPECTROGRAM_HPP
