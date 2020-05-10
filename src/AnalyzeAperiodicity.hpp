// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_DSP_ANALYZE_APERIODICITY_HPP
#define UZUME_DSP_ANALYZE_APERIODICITY_HPP
#include "Spectrum.hpp"

namespace uzume { namespace dsp {

class InstantWaveform {
public:
    const double *wave;
    const double f0;
    const unsigned int length;
    const unsigned int samplingFrequency;
};

class AnalyzeAperiodicity {
public:
    virtual ~AnalyzeAperiodicity() = default;

    virtual bool operator()(Spectrum *output, const InstantWaveform *input) = 0;
};

} }

#endif //UZUME_DSP_ANALYZE_APERIODICITY_HPP
