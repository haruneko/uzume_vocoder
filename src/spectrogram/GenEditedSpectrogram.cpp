// Copyright 2022 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include <cmath>
#include "GenEditedSpectrogram.hpp"

using namespace uzume::vocoder;

GenEditedSpectrogram::GenEditedSpectrogram(const Spectrogram *spectrogram, const ControlChange &cc)
        : spectrogram(spectrogram), cc(cc) {
}

bool GenEditedSpectrogram::pickUpSpectrumAt(Spectrum *destination, double ms) const {
    Spectrum s(destination->fftSize);
    bool result = spectrogram->pickUpSpectrumAt(&s, ms);
    if(!result) return false;

    double ratio = pow(2, cc.at(ms / msLength()));
    int maxIndex = fftSize() / 2;
    for(int i = 0; i <= destination->fftSize / 2; i++) {
        double dIndex = (double)i * ratio;
        int indexFloor = (int)floor(dIndex);
        int indexCeil = indexFloor + 1;
        double interop = dIndex - indexFloor;
        destination->periodicSpectrum[i] =
                s.periodicSpectrum[std::min(maxIndex, indexFloor)] * (1 - interop) +
                s.periodicSpectrum[std::min(maxIndex, indexCeil)] * interop;
        destination->aperiodicSpectrum[i] =
                s.aperiodicSpectrum[std::min(maxIndex, indexFloor)] * (1 - interop) +
                s.aperiodicSpectrum[std::min(maxIndex, indexCeil)] * interop;
    }
    return result;
}

double GenEditedSpectrogram::f0At(double ms) const {
    return spectrogram->f0At(ms);
}

double GenEditedSpectrogram::msLength() const {
    return spectrogram->msLength();
}

unsigned int GenEditedSpectrogram::fftSize() const {
    return spectrogram->fftSize();
}
