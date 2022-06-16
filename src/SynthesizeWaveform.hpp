// Copyright 2021 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef __UZUME_VOCODER_SYNTHESIZE_WAVEFORM_HPP__
#define __UZUME_VOCODER_SYNTHESIZE_WAVEFORM_HPP__
#include "Spectrogram.hpp"
#include "data/Waveform.hpp"

namespace uzume { namespace vocoder {

/**
 * SynthesizeWaveform is an interface that synthesize waveform from given spectrogram.
 */
class SynthesizeWaveform {
public:
    SynthesizeWaveform() = default;
    SynthesizeWaveform(const SynthesizeWaveform&) = delete;
    virtual ~SynthesizeWaveform() = default;

    /**
     * operator()
     * @param output is a target waveform to write synthesized waveform.
     * @param input is a source to synthesize from.
     * @return true when synthesis succeeds, otherwise false.
     */
    virtual bool operator()(Waveform *output, const Spectrogram *input) = 0;
};

} }

#endif // __UZUME_VOCODER_SYNTHESIZE_WAVEFORM_HPP__
