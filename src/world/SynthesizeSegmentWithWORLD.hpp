// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef __UZUME_VOCODER_SYNTHESIZE_SEGMENT_WITH_WORLD_HPP__
#define __UZUME_VOCODER_SYNTHESIZE_SEGMENT_WITH_WORLD_HPP__
#include "../SynthesizeImpulseResponse.hpp"
#include "../SynthesizeSegment.hpp"

namespace uzume { namespace vocoder { namespace world {

/**
 * SynthesizeSegmentWithWORLD is an implementation of SynthesizeSegment with WORLD.
 */
class SynthesizeSegmentWithWORLD final : public SynthesizeSegment {
public:
    SynthesizeSegmentWithWORLD() = delete;
    explicit SynthesizeSegmentWithWORLD(SynthesizeImpulseResponse *synthesize);
    ~SynthesizeSegmentWithWORLD() override = default;

    bool operator()(SegmentSignal *output, SegmentParameters *input) override;
private:
    SynthesizeImpulseResponse *synthesize;
};

} } }

#endif // __UZUME_VOCODER_SYNTHESIZE_SEGMENT_WITH_WORLD_HPP__
