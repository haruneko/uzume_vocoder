// Copyright 2021 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef __UZUME_VOCODER_SYNTHESIZE_WAVEFORM_WITH_WORLD_HPP__
#define __UZUME_VOCODER_SYNTHESIZE_WAVEFORM_WITH_WORLD_HPP__

#include <functional>

#include "../SynthesizeWaveform.hpp"

namespace uzume { namespace vocoder {
class SynthesizeImpulseResponse;
class SynthesizeSegment;
} }

namespace uzume { namespace vocoder { namespace world {

/**
 * SynthesizeWaveformWithWORLD is an implementation of SynthesizeWaveform with WORLD.
 */
class SynthesizeWaveformWithWORLD final : public SynthesizeWaveform {
public:
    static const std::function<SynthesizeImpulseResponse *(unsigned int, unsigned int)> DefaultSynthesizeImpulseResponseFactory;
    static const std::function<SynthesizeSegment *(SynthesizeImpulseResponse *)> DefaultSynthesizeSegmentFactory;

    explicit SynthesizeWaveformWithWORLD(
        const std::function<SynthesizeImpulseResponse *(unsigned int, unsigned int)> &synthesizeImpulseResponseFactory = DefaultSynthesizeImpulseResponseFactory,
        const std::function<SynthesizeSegment *(SynthesizeImpulseResponse *)> &synthesizeSegmentFactory = DefaultSynthesizeSegmentFactory);
    bool operator()(Waveform *output, const Spectrogram *input) override;

private:
    const std::function<SynthesizeImpulseResponse *(unsigned int, unsigned int)> synthesizeImpulseResponseFactory;
    const std::function<SynthesizeSegment *(SynthesizeImpulseResponse *)> synthesizeSegmentFactory;
};

} } }

#endif // __UZUME_VOCODER_SYNTHESIZE_WAVEFORM_WITH_WORLD_HPP__
