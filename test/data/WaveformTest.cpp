// Copyright 2021 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "data/Waveform.hpp"
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class WaveformTest : public ::testing::Test {
public:

    virtual void SetUp() override {
        waveform = new uzume::vocoder::Waveform(44100, 44100);
        // initialize with 440 Hz sin wave.
        for(int i = 0; i < 44100; i++) {
            waveform->data[i] = sin(2 * M_PI * (double)i / 44100 * 440.0);
        }
    }

    virtual void TearDown() override {
        delete waveform;
    }

    uzume::vocoder::Waveform *waveform;
};

using namespace uzume::vocoder;

TEST_F(WaveformTest, constructorShouldClearAllDataToZero) {
    auto *waveform = new uzume::vocoder::Waveform(44100, 44100);
    for(unsigned int i = 0; i < waveform->length; i++) {
        EXPECT_DOUBLE_EQ(waveform->data[i], 0.0);
    }
    delete waveform;
}

TEST_F(WaveformTest, msLengthShouldReturnMilliSeconds) {
    EXPECT_DOUBLE_EQ(waveform->msLength(), 1000.0);
}

TEST_F(WaveformTest, indexAtShouldReturnTheCorrespondingIndex) {
    EXPECT_EQ(waveform->indexAt(-500.0), -22050);
    EXPECT_EQ(waveform->indexAt(0.0), 0);
    EXPECT_EQ(waveform->indexAt(500.0), 22050);
    EXPECT_EQ(waveform->indexAt(1000.0), 44100);
    EXPECT_EQ(waveform->indexAt(2000.0), 88200);
}

TEST_F(WaveformTest, maxAbsoluteValueBetweenShouldReturnTheMaxAbsoluteValue) {
    EXPECT_DOUBLE_EQ(0.0, waveform->maxAbsoluteValueBetween(-1000.0, 0.0));
    EXPECT_DOUBLE_EQ(0.0, waveform->maxAbsoluteValueBetween(0.0, 0.0));
    EXPECT_DOUBLE_EQ(0.0, waveform->maxAbsoluteValueBetween(1000.0, 2000.0));
    EXPECT_NEAR(1.0, waveform->maxAbsoluteValueBetween(0.0, 1000.0), 0.00001);
}

TEST_F(WaveformTest, rootMeanSquareBetweenShouldReturnReturnRMS) {
    EXPECT_DOUBLE_EQ(0.0, waveform->rootMeanSquareBetween(-1000.0, 0.0));
    EXPECT_DOUBLE_EQ(0.0, waveform->rootMeanSquareBetween(0.0, 0.0));
    EXPECT_DOUBLE_EQ(0.0, waveform->rootMeanSquareBetween(1000.0, 2000.0));
    EXPECT_DOUBLE_EQ(0.70710678118654757, waveform->rootMeanSquareBetween(0.0, 1000.0));
}