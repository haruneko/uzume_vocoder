// Copyright 2022 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_VOCODER_GEN_EDITED_SPECTROGRAM_HPP
#define UZUME_VOCODER_GEN_EDITED_SPECTROGRAM_HPP

#include "../Spectrogram.hpp"
#include "../data/ControlChange.hpp"

namespace uzume { namespace vocoder {

class GenEditedSpectrogram final : public Spectrogram {
public:
    GenEditedSpectrogram(const Spectrogram *spectrogram, const ControlChange &cc);

    bool pickUpSpectrumAt(Spectrum *destination, double ms) const override;
    double f0At(double ms) const override;
    double msLength() const override;
    unsigned int fftSize() const override;
private:
    const Spectrogram *spectrogram;
    ControlChange cc;
};

} }

#endif //UZUME_VOCODER_GEN_EDITED_SPECTROGRAM_HPP