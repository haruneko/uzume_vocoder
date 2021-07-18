// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_VOCODER_WAVEFORM_SPECTROGRAM_HPP
#define UZUME_VOCODER_WAVEFORM_SPECTROGRAM_HPP

#include <functional>

#include "../AnalyzeAperiodicity.hpp"
#include "../AnalyzePeriodicity.hpp"
#include "../data/Contour.hpp"
#include "../EstimateF0.hpp"
#include "../data/InstantWaveform.hpp"
#include "../Spectrogram.hpp"
#include "../data/Waveform.hpp"

namespace uzume { namespace vocoder {

/**
 * WaveformSpectrogram is an implementation of Spectrogram.
 * This class treat waveform as spectrogram with WORLD analysis algorithm.
 */
class WaveformSpectrogram final : public Spectrogram {
public:
    static const std::function<AnalyzeAperiodicity *(unsigned int)> DefaultAperiodicAnalysisFactory;
    static const std::function<AnalyzePeriodicity *(unsigned int)> DefaultPeriodicAnalysisFactory;
    static const std::function<EstimateF0 *(double)> DefaultF0EstimationFactory;

    explicit WaveformSpectrogram(Waveform *waveform,
                                 const std::function<AnalyzeAperiodicity *(
                                         unsigned int)> &aperiodicAnalysisFactory = DefaultAperiodicAnalysisFactory,
                                 const std::function<AnalyzePeriodicity *(
                                         unsigned int)> &periodicAnalysisFactory = DefaultPeriodicAnalysisFactory,
                                 const std::function<EstimateF0 *(
                                         double)> &f0EstimationFactory = DefaultF0EstimationFactory);

    ~WaveformSpectrogram() final;

    bool pickUpSpectrumAt(Spectrum *destination, double ms) const override;

    double f0At(double ms) const override;

    double msLength() const override;

    unsigned int fftSize() const override;

private:
    Waveform *waveform; // Note that waveform is `not` generated in this class.

    AnalyzeAperiodicity *analyzeAperiodicity;
    AnalyzePeriodicity *analyzePeriodicity;
    Contour *f0;
    InstantWaveform *iw;
};

} }

#endif //UZUME_VOCODER_WAVEFORM_SPECTROGRAM_HPP
