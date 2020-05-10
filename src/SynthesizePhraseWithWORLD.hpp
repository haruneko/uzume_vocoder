// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_DSP_SYNTHESIZE_PHRASE_WITH_WORLD_HPP
#define UZUME_DSP_SYNTHESIZE_PHRASE_WITH_WORLD_HPP
#include "SynthesizeImpulseResponse.hpp"
#include "SynthesizePhrase.hpp"

namespace uzume { namespace dsp {

class SynthesizePhraseWithWORLD final : public SynthesizePhrase {
public:
    SynthesizePhraseWithWORLD() = delete;
    explicit SynthesizePhraseWithWORLD(SynthesizeImpulseResponse *synthesize);
    ~SynthesizePhraseWithWORLD() override = default;

    bool operator()(PhraseSignal *output, PhraseParameters *input) override;
private:
    SynthesizeImpulseResponse *synthesize;
};

} }

#endif //UZUME_DSP_SYNTHESIZE_PHRASE_WITH_WORLD_HPP
