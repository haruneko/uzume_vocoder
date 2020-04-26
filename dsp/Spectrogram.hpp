// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef __DSP_SPECTROGRAM_HPP__
#define __DSP_SPECTROGRAM_HPP__

#include "Spectrum.hpp"

/**
 * Spectrogram represents an interface of vocal spectrogram.
 */
class Spectrogram {
public:
    virtual ~Spectrogram() = default;

    /**
     * Spectrogram#pickUpSpectrumAt picks up an instant spectrum on the specified moment.
     * @param destination is an output buffer of spectrum.
     * @param ms is the time to pick up.
     * @return true if picking up is successful.
     * @return false otherwise.
     */
    virtual bool pickUpSpectrumAt(Spectrum *destination, double ms) const = 0;

    /**
     * Spectrogram#f0At returns a f0 value of spectrogram at ms[milli seconds].
     */
    virtual double f0At(double ms) const = 0;

    /**
     * Spectrogram#msLength returns a length of spectrogram in milli seconds.
     */
    virtual double msLength() const = 0;

    /**
     * Spectrogram#fftLength returns a FFT length of spectrogram.
     */
    virtual unsigned int fftSize() const = 0;
};

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

#endif

