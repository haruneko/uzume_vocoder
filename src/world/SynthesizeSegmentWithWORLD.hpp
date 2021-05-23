// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_DSP_SYNTHESIZE_Segment_WITH_WORLD_HPP
#define UZUME_DSP_SYNTHESIZE_Segment_WITH_WORLD_HPP
#include "../SynthesizeImpulseResponse.hpp"
#include "../SynthesizeSegment.hpp"

namespace uzume { namespace dsp { namespace world {

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

#endif //UZUME_DSP_SYNTHESIZE_Segment_WITH_WORLD_HPP
