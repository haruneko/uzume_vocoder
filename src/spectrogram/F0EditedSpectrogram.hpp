// Copyright 2022 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_VOCODER_FO_EDITED_SPECTROGRAM_HPP
#define UZUME_VOCODER_FO_EDITED_SPECTROGRAM_HPP

#include <functional>
#include "../Spectrogram.hpp"
#include "../data/ControlChange.hpp"

namespace uzume { namespace vocoder {

enum class SynthType {
    Linear,
    Log
};


/**
 * F0EditedSpectrogram is a wrapper spectrogram with F0 edited by ControlChange.
 * Everything except f0 is all the same as the given spectrogram,
 * but f0 is multiplied by the value, ratio in the given ControlChange.
 */
class F0EditedSpectrogram final : public Spectrogram {
public:
    F0EditedSpectrogram(const Spectrogram *spectrogram, const ControlChange &cc, const SynthType synthType);

    bool pickUpSpectrumAt(Spectrum *destination, double ms) const override;
    double f0At(double ms) const override;
    double msLength() const override;
    unsigned int fftSize() const override;
private:
    const Spectrogram *spectrogram;
    ControlChange cc;
    std::function<double(double, double)> f;
};

} }


#endif //UZUME_VOCODER_FO_EDITED_SPECTROGRAM_HPP