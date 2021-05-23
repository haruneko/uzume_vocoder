// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef __DSP_SYNTHESIZE_Segment_HPP__
#define __DSP_SYNTHESIZE_Segment_HPP__

#include "Spectrogram.hpp"

namespace uzume { namespace dsp {

/**
 * SegmentSignal is a value class for Segment synthesis output.
 * This class is not responsible for memory allocation.
 * Make sure neither SegmentSignal#indexMin nor SegmentSignal#indexMax is out of array range.
 */
class SegmentSignal final {
public:
    SegmentSignal() = delete;
    SegmentSignal(double *waveStartPoint, int indexMin, int indexMax, unsigned int samplingFrequency);
    double *waveStartPoint;
    int indexMin;
    int indexMax;
    unsigned int samplingFrequency;
};

/**
 * SegmentSignal is a value class for Segment synthesis input.
 * This class is not responsible for memory allocation.
 * SegmentParameters#startPhase and SegmentParameters#startFractionalTimeShift should be set
 * when SynthesizeSegment begins synthesis from the middle of music sentence, otherwise they should be 0.0.
 */
class SegmentParameters final {
public:
    SegmentParameters() = delete;
    SegmentParameters(Spectrogram *spectrogram, double startPhase, double startFractionalTimeShift);

    Spectrogram *spectrogram;
    double startPhase;
    double startFractionalTimeShift;
};

/**
 * SynthesizeSegment synthesize a single Segment.
 */
class SynthesizeSegment {
public:
    SynthesizeSegment() = default;
    SynthesizeSegment(const SynthesizeSegment&) = delete;
    virtual ~SynthesizeSegment() = default;

    /**
     * SynthesizeSegment#() synthesize a single Segment into output from input.
     * @return true if synthesis is successful.
     * @return false otherwise.
     */
    virtual bool operator()(SegmentSignal *output, SegmentParameters *input) = 0;
};

} }

#endif
