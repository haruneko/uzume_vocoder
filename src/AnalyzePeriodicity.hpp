// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_DSP_ANALYZE_PERIODICITY_HPP
#define UZUME_DSP_ANALYZE_PERIODICITY_HPP
#include "data/InstantWaveform.hpp"
#include "data/Spectrum.hpp"

namespace uzume { namespace dsp {

/**
 * AnalyzePeriodicity is an interface to analyze periodic spectrum.
 */
class AnalyzePeriodicity {
public:
    virtual ~AnalyzePeriodicity() = default;

    /**
     * () analyzes instant waveform in input and sets periodic spectrum into output.
     * `input` should be instantaneous signal that contains, at least, a single period of voice.
     * @return true if analysis is successful, otherwise false.
     */
    virtual bool operator()(Spectrum *output, const InstantWaveform *input) = 0;
};

} }

#endif //UZUME_DSP_ANALYZE_PERIODICITY_HPP
