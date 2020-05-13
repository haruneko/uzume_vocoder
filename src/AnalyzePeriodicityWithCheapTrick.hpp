// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_DSP_ANALYZE_PERIODICITY_WITH_CHEAPTRICK_HPP
#define UZUME_DSP_ANALYZE_PERIODICITY_WITH_CHEAPTRICK_HPP
#include "fft.hpp"
#include "AnalyzePeriodicity.hpp"
#include "GaussianNoiseGenerator.hpp"

namespace uzume { namespace dsp {

class AnalyzePeriodicityWithCheapTrick final : public AnalyzePeriodicity {
public:
    AnalyzePeriodicityWithCheapTrick() = delete;
    AnalyzePeriodicityWithCheapTrick(unsigned int fftSize, unsigned int samplingFrequency);
    ~AnalyzePeriodicityWithCheapTrick();

    bool operator()(Spectrum *output, const InstantWaveform *input) override;

    const unsigned int fftSize;
    const unsigned int samplingFrequency;

private:
    ForwardRealFFT forwardRealFft;
    InverseRealFFT inverseRealFft;

    GaussianNoiseGenerator randn;
};

} }

#endif //UZUME_DSP_ANALYZE_PERIODICITY_WITH_CHEAPTRICK_HPP
