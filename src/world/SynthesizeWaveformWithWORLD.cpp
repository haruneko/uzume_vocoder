// Copyright 2021 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include "../Spectrogram.hpp"
#include "../data/Waveform.hpp"

#include "SynthesizeImpulseResponseWithWORLD.hpp"
#include "SynthesizeSegmentWithWORLD.hpp"
#include "SynthesizeWaveformWithWORLD.hpp"

using namespace uzume::vocoder;
using namespace uzume::vocoder::world;

const std::function<SynthesizeImpulseResponse *(unsigned int, unsigned int)>
        SynthesizeWaveformWithWORLD::DefaultSynthesizeImpulseResponseFactory = [](unsigned int fftSize, unsigned int samplingFrequency) {
    return new SynthesizeImpulseResponseWithWORLD(fftSize, samplingFrequency);
};

const std::function<SynthesizeSegment *(SynthesizeImpulseResponse *)>
        SynthesizeWaveformWithWORLD::DefaultSynthesizeSegmentFactory = [](SynthesizeImpulseResponse *irs) {
    return new SynthesizeSegmentWithWORLD(irs);
};

SynthesizeWaveformWithWORLD::SynthesizeWaveformWithWORLD(
        const std::function<SynthesizeImpulseResponse *(unsigned int, unsigned int)> &synthesizeImpulseResponseFactory,
        const std::function<SynthesizeSegment *(SynthesizeImpulseResponse *)> &synthesizeSegmentFactory) :
        synthesizeImpulseResponseFactory(synthesizeImpulseResponseFactory),
        synthesizeSegmentFactory(synthesizeSegmentFactory) {
}

bool SynthesizeWaveformWithWORLD::operator()(Waveform *output, const Spectrogram *input) {
    if(!output || output->length == 0 || !input || input->msLength() <= 0.0) {
        return false;
    }

    auto irs = synthesizeImpulseResponseFactory(input->fftSize(), output->samplingFrequency);
    auto ss = synthesizeSegmentFactory(irs);

    for(unsigned int i = 0; i < output->length; i++) {
        output->data[i] = 0.0;
    }

    SegmentSignal s(output->data, /* indexMin = */ 0, /* indexMax = */ output->length, output->samplingFrequency);
    SegmentParameters p(input, /* startPhase = */ 0.0, /* startFractionalTimeShift = */ 0.0);

    (*ss)(&s, &p);

    delete ss;
    delete irs;

    return true;
}
