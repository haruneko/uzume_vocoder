// Copyright 2021 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "../MockSpectrogram.hpp"
#include "../MockTimeAxisMap.hpp"
#include "spectrogram/StretchedPartialSpectrogram.hpp"

class StretchedPartialSpectrogramTest : public ::testing::Test {
public:
    virtual void setUp() { }
};

using namespace uzume::vocoder;

TEST_F(StretchedPartialSpectrogramTest, pickUpSpectrumAtShouldReturnSpectrumAtTheTimeTimeAxisMapPoints) {
    Spectrum s(1024);
    MockTimeAxisMap tam;
    EXPECT_CALL(tam, at(2.0)).Times(1).WillOnce(::testing::Return(100.0));
    MockSpectrogram sp;
    EXPECT_CALL(sp, pickUpSpectrumAt(&s, 100.0)).Times(1).WillOnce(::testing::Return(true));

    StretchedPartialSpectrogram sut(&sp, &tam);

    EXPECT_TRUE(sut.pickUpSpectrumAt(&s, 2.0));
}

TEST_F(StretchedPartialSpectrogramTest, f0AtShouldReturnF0AtTheTimeTimeAxisMapPoints) {
    MockTimeAxisMap tam;
    EXPECT_CALL(tam, at(2.0)).Times(1).WillOnce(::testing::Return(100.0));
    MockSpectrogram sp;
    EXPECT_CALL(sp, f0At(100.0)).Times(1).WillOnce(::testing::Return(440.0));

    StretchedPartialSpectrogram sut(&sp, &tam);

    EXPECT_DOUBLE_EQ(sut.f0At(2.0), 440.0);
}

TEST_F(StretchedPartialSpectrogramTest, msLengthShouldUseTimeAxisValue) {
    MockTimeAxisMap tam;
    EXPECT_CALL(tam, msLength()).Times(1).WillOnce(::testing::Return(1234.0));

    StretchedPartialSpectrogram sut(nullptr, &tam);

    EXPECT_DOUBLE_EQ(sut.msLength(), 1234.0);
}

TEST_F(StretchedPartialSpectrogramTest, fftSizeShouldDelegateDependentSpectrogramsFftSize) {
    MockSpectrogram sp;
    EXPECT_CALL(sp, fftSize()).Times(1).WillOnce(::testing::Return(12345));

    StretchedPartialSpectrogram sut(&sp, nullptr);

    EXPECT_EQ(sut.fftSize(), 12345);
}
