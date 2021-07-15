// Copyright 2021 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "../MockSpectrogram.hpp"
#include "spectrogram/ArraySpectrogramAggregator.hpp"

class ArraySpectrogramAggregatorTest : public ::testing::Test {
public:
    virtual void setUp() { }
};

using namespace uzume::vocoder;

TEST_F(ArraySpectrogramAggregatorTest, msLengthShouldReturnSummaryLengthOfAggregatedSpectrogram) {
    MockSpectrogram sp1, sp2, sp3;
    double sp1Length = 0.125, sp2Length = 0.5, sp3Length = 1.0;
    EXPECT_CALL(sp1, msLength()).Times(1).WillOnce(::testing::Return(sp1Length));
    EXPECT_CALL(sp2, msLength()).Times(1).WillOnce(::testing::Return(sp2Length));
    EXPECT_CALL(sp3, msLength()).Times(1).WillOnce(::testing::Return(sp3Length));
    EXPECT_CALL(sp1, fftSize()).WillRepeatedly(::testing::Return(1024));
    EXPECT_CALL(sp2, fftSize()).WillRepeatedly(::testing::Return(1024));
    EXPECT_CALL(sp3, fftSize()).WillRepeatedly(::testing::Return(1024));

    std::vector<const Spectrogram *> a = {&sp1, &sp2, &sp3};
    auto s = ArraySpectrogramAggregator::from(a);

    EXPECT_DOUBLE_EQ(s->msLength(), sp1Length + sp2Length + sp3Length);
    delete s;
}

TEST_F(ArraySpectrogramAggregatorTest, f0AtShouldReturnTheCorrespondingF0InArray) {
    MockSpectrogram sp1, sp2, sp3;
    double sp1Length = 100.0, sp2Length = 50.0, sp3Length = 50.0;
    EXPECT_CALL(sp1, msLength()).WillRepeatedly(::testing::Return(sp1Length));
    EXPECT_CALL(sp2, msLength()).WillRepeatedly(::testing::Return(sp2Length));
    EXPECT_CALL(sp3, msLength()).WillRepeatedly(::testing::Return(sp3Length));
    EXPECT_CALL(sp1, fftSize()).WillRepeatedly(::testing::Return(1024));
    EXPECT_CALL(sp2, fftSize()).WillRepeatedly(::testing::Return(1024));
    EXPECT_CALL(sp3, fftSize()).WillRepeatedly(::testing::Return(1024));
    EXPECT_CALL(sp2, f0At(125.0 - 100.0)).WillOnce(::testing::Return(220.0));
    EXPECT_CALL(sp3, f0At(160.0 - 150.0)).WillOnce(::testing::Return(440.0));

    std::vector<const Spectrogram *> a = {&sp1, &sp2, &sp3};
    auto s = ArraySpectrogramAggregator::from(a);

    EXPECT_DOUBLE_EQ(s->f0At(125.0), 220.0);
    EXPECT_DOUBLE_EQ(s->f0At(160.0), 440.0);
    delete s;
}

TEST_F(ArraySpectrogramAggregatorTest, pickUpSpectrumAtShouldPickUpSpectrumAtCorrespondingTime) {
    MockSpectrogram sp1, sp2, sp3;
    Spectrum ss(1024);
    double sp1Length = 100.0, sp2Length = 50.0, sp3Length = 50.0;
    EXPECT_CALL(sp1, msLength()).WillRepeatedly(::testing::Return(sp1Length));
    EXPECT_CALL(sp2, msLength()).WillRepeatedly(::testing::Return(sp2Length));
    EXPECT_CALL(sp3, msLength()).WillRepeatedly(::testing::Return(sp3Length));
    EXPECT_CALL(sp1, fftSize()).WillRepeatedly(::testing::Return(1024));
    EXPECT_CALL(sp2, fftSize()).WillRepeatedly(::testing::Return(1024));
    EXPECT_CALL(sp3, fftSize()).WillRepeatedly(::testing::Return(1024));
    EXPECT_CALL(sp2, pickUpSpectrumAt(&ss, 125.0 - 100.0)).WillOnce(::testing::Return(true));
    EXPECT_CALL(sp3, pickUpSpectrumAt(&ss, 160.0 - 150.0)).WillOnce(::testing::Return(false));

    std::vector<const Spectrogram *> a = {&sp1, &sp2, &sp3};
    auto s = ArraySpectrogramAggregator::from(a);

    EXPECT_EQ(s->pickUpSpectrumAt(&ss, 125.0), true);
    EXPECT_EQ(s->pickUpSpectrumAt(&ss, 160.0), false);
    delete s;
}
