// Copyright 2021 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_VOCODER_STRETCHED_PARTIAL_SPECTROGRAM_HPP
#define UZUME_VOCODER_STRETCHED_PARTIAL_SPECTROGRAM_HPP

#include "../Spectrogram.hpp"
#include "../TimeAxisMap.hpp"

namespace uzume { namespace vocoder {

class StretchedPartialSpectrogram final : public Spectrogram {
public:
    StretchedPartialSpectrogram() = delete;
    StretchedPartialSpectrogram(const Spectrogram *spectrogram, const TimeAxisMap *timeAxisMap);
    ~StretchedPartialSpectrogram() final = default;

    bool pickUpSpectrumAt(Spectrum *destination, double ms) const final;
    double f0At(double ms) const final;
    double msLength() const final;
    unsigned int fftSize() const final;

private:
    const Spectrogram *spectrogram;
    const TimeAxisMap *timeAxisMap;
};

} }

#endif // UZUME_VOCODER_STRETCHED_PARTIAL_SPECTROGRAM_HPP