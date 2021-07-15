// Copyright 2021 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef __UZUME_VOCODER_ARRAY_SPECTROGRAM_AGGREGATOR__
#define __UZUME_VOCODER_ARRAY_SPECTROGRAM_AGGREGATOR__
#include <vector>

#include "../Spectrogram.hpp"

namespace uzume { namespace vocoder {

/**
 * ArraySpectrogramAggregator is an aggregation of Spectrograms.
 * ArraySpectrogramAggregator simply connects all spectrograms in series.
 */
class ArraySpectrogramAggregator final : public Spectrogram {
public:
    ArraySpectrogramAggregator() = delete;

    static ArraySpectrogramAggregator *from(std::vector<const Spectrogram *> &spectrograms);

    double f0At(double ms) const override;
    bool pickUpSpectrumAt(Spectrum *destination, double ms) const override;
    double msLength() const override;
    unsigned int fftSize() const override;

private:
    explicit ArraySpectrogramAggregator(std::vector<const Spectrogram *> &spectrograms);

    /**
     * @param ms represents a position in a entire Spectrogram that connects each element in array.
     * @return corresponding index and ms in partial spectrogram at `ms`.
     */
    std::pair<int, double> indexAndMsAt(double ms) const;

    std::vector<const Spectrogram *> spectrograms;
};
}}

#endif // __UZUME_VOCODER_ARRAY_SPECTROGRAM_AGGREGATOR__
