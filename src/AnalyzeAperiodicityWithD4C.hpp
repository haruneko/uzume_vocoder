// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_DSP_ANALYZE_APERIODICITY_WITH_D4C_HPP
#define UZUME_DSP_ANALYZE_APERIODICITY_WITH_D4C_HPP

#include "AnalyzeAperiodicity.hpp"
#include "fft.hpp"
#include "GaussianNoiseGenerator.hpp"
#include "Spectrum.hpp"

namespace uzume { namespace dsp {

/**
 * AnalyzeAperiodicityWithD4C is an implementation of aperiodic spectrum analysis.
 */
class AnalyzeAperiodicityWithD4C final : public AnalyzeAperiodicity {
public:
    AnalyzeAperiodicityWithD4C() = delete;
    AnalyzeAperiodicityWithD4C(unsigned int fftSize, unsigned int samplingFrequency);
    ~AnalyzeAperiodicityWithD4C();

    /**
     * () analyzes input with D4C and sets its aperiodic spectrum into output.
     */
    bool operator()(Spectrum *output, const InstantWaveform *input) override;

    const unsigned int fftSize;
    const unsigned int samplingFrequency;

private:
    int numberOfAperiodicities() const;
    int fftSizeForD4C() const;
    int nuttallWindowSize() const;
    int fftSizeForD4CLoveTrain() const;

    ForwardRealFFT forwardRealFft;
    InverseRealFFT inverseRealFft;
    double *coarseAperiodicity;
    double *coarseFrequencyAxis;
    double *frequencyAxis;
    double *nuttallWindow;

    ForwardRealFFT forwardRealFFtForD4CLoveTrain;
    const int boundary0;
    const int boundary1;
    const int boundary2;

    GaussianNoiseGenerator randn;
};

} }

#endif //UZUME_DSP_ANALYZE_APERIODICITY_WITH_D4C_HPP
