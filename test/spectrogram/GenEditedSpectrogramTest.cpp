// Copyright 2022 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "../MockSpectrogram.hpp"
#include "spectrogram/GenEditedSpectrogram.hpp"

using namespace uzume::vocoder;

class GenEditedSpectrogramTest : public ::testing::Test {
public:
    virtual void setUp() { }
};

TEST_F(GenEditedSpectrogramTest, pickUpSpectrumAtShouldStretchPickedUpSpectrum) {
    MockSpectrogram sp;
    Spectrum ss(1024);

    EXPECT_CALL(sp, msLength()).WillRepeatedly(::testing::Return(500.0));
    EXPECT_CALL(sp, fftSize()).WillRepeatedly(::testing::Return(1024));
    EXPECT_CALL(sp, pickUpSpectrumAt(::testing::_, ::testing::_)).WillRepeatedly(::testing::Return(true));

    auto s = new GenEditedSpectrogram(&sp, ControlChange(1.25));

    EXPECT_EQ(s->pickUpSpectrumAt(&ss, 125.0), true);
    delete s;
}
