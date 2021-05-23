// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include <gtest/gtest.h>

#include "data/Waveform.hpp"
#include "spectrogram/WaveformSpectrogram.hpp"

class WaveformSpectrogramTest : public ::testing::Test {
protected:
    virtual void setUp() { }
};

TEST_F(WaveformSpectrogramTest, pickUpSpectrumAtShouldNotThrowSIGSEGV) {
    auto *w = new uzume::dsp::Waveform(44100, 44100);
    uzume::dsp::WaveformSpectrogram ws(w);
    uzume::dsp::Spectrum s(ws.fftSize());

    for(int i = -5; i < 15; i++) {
        EXPECT_TRUE(ws.pickUpSpectrumAt(&s, i * 100.0));
    }
}
