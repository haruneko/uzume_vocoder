// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include "SynthesizeSegment.hpp"

using namespace uzume::vocoder;

SegmentSignal::SegmentSignal(double *waveStartPoint, int indexMin, int indexMax, unsigned int samplingFrequency)
        : waveStartPoint(waveStartPoint), indexMin(indexMin), indexMax(indexMax), samplingFrequency(samplingFrequency) {
}

SegmentParameters::SegmentParameters(Spectrogram *spectrogram, double startPhase, double startFractionalTimeShift)
        : spectrogram(spectrogram), startPhase(startPhase), startFractionalTimeShift(startFractionalTimeShift) {
}
