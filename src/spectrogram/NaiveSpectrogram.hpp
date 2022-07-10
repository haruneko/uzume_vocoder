// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_VOCODER_NAIVE_SPECTROGRAM_HPP
#define UZUME_VOCODER_NAIVE_SPECTROGRAM_HPP

#include "../Spectrogram.hpp"
#include "../data/Contour.hpp"
#include "../data/Waveform.hpp"

namespace uzume { namespace vocoder {

/**
 * NaiveSpectrogram represents a Spectrogram that
 * contains raw f0Contour, PeriodicSpectrogram and AperiodicSpectrogram.
 */
class NaiveSpectrogram final : public Spectrogram {
public:
    NaiveSpectrogram() = delete;
    NaiveSpectrogram(unsigned int length, unsigned int fftSize, double msFramePeriod);

    ~NaiveSpectrogram() override;

    bool pickUpSpectrumAt(Spectrum *destination, double ms) const override final;
    double f0At(double ms) const override final;
    double msLength() const override final;
    unsigned int fftSize() const override final;

    static NaiveSpectrogram *from(const Waveform *wave, double msFramePeriod);

    double **periodicSpecgram;
    double **aperiodicSpecgram;
    Contour *f0Contour;
    double *timeAxis;

private:
    unsigned int length;
    unsigned int _fftSize;
    double msFramePeriod;
};

} }

#endif //UZUME_VOCODER_NAIVE_SPECTROGRAM_HPP
