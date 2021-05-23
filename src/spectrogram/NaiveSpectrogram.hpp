// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_DSP_NAIVE_SPECTROGRAM_HPP
#define UZUME_DSP_NAIVE_SPECTROGRAM_HPP
#include "../Spectrogram.hpp"

namespace uzume { namespace dsp {

class NaiveSpectrogram final : public Spectrogram {
public:
    NaiveSpectrogram() = delete;
    NaiveSpectrogram(unsigned int length, unsigned int fftSize, double msFramePeriod);

    ~NaiveSpectrogram() override;

    bool pickUpSpectrumAt(Spectrum *destination, double ms) const override;
    double f0At(double ms) const override;
    double msLength() const override;
    unsigned int fftSize() const override;

    double **periodicSpecgram;
    double **aperiodicSpecgram;
    double *f0Contour;
    double *timeAxis;

private:
    unsigned int length;
    unsigned int _fftSize;
    double msFramePeriod;
};

} }

#endif //UZUME_DSP_NAIVE_SPECTROGRAM_HPP
