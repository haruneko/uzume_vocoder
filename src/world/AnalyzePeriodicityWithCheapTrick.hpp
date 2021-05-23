// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_VOCODER_ANALYZE_PERIODICITY_WITH_CHEAPTRICK_HPP
#define UZUME_VOCODER_ANALYZE_PERIODICITY_WITH_CHEAPTRICK_HPP
#include "fft.hpp"
#include "../AnalyzePeriodicity.hpp"
#include "GaussianNoiseGenerator.hpp"

namespace uzume { namespace vocoder { namespace world {

/**
 * AnalyzePeriodicityWithCheapTrick is an implementation of periodic spectrum analysis.
 */
class AnalyzePeriodicityWithCheapTrick final : public AnalyzePeriodicity {
public:
    AnalyzePeriodicityWithCheapTrick() = delete;
    explicit AnalyzePeriodicityWithCheapTrick(unsigned int samplingFrequency);
    ~AnalyzePeriodicityWithCheapTrick() final;

    /**
     * () analyze input with Cheap Trick and puts its periodic spectrum into output.
     */
    bool operator()(Spectrum *output, const InstantWaveform *input) override;

    const unsigned int fftSize;
    const unsigned int samplingFrequency;

private:
    ForwardRealFFT forwardRealFft;
    InverseRealFFT inverseRealFft;

    GaussianNoiseGenerator randn;
};

} } }

#endif //UZUME_VOCODER_ANALYZE_PERIODICITY_WITH_CHEAPTRICK_HPP
