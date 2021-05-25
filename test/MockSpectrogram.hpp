// Copyright 2021 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_VOCODER_MOCK_SPECTROGRAM_HPP
#define UZUME_VOCODER_MOCK_SPECTROGRAM_HPP

#include <gmock/gmock.h>
#include "Spectrogram.hpp"

namespace uzume { namespace vocoder {
class MockSpectrogram : public Spectrogram {
public:
    MOCK_METHOD(bool, pickUpSpectrumAt, (Spectrum *, double ), (const, override));
    MOCK_METHOD(double, f0At, (double), (const, override));
    MOCK_METHOD(double, msLength, (), (const, override));
    MOCK_METHOD(unsigned int, fftSize, (), (const, override));
};
} }

#endif // UZUME_VOCODER_MOCK_SPECTROGRAM_HPP
