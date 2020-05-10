// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_DSP_ESTIMATE_F0_HPP
#define UZUME_DSP_ESTIMATE_F0_HPP
#include "Waveform.hpp"

namespace uzume { namespace dsp {

class F0Contour final {
public:
    double *f0;
    int length;
    double msFramePeriod;
};

class EstimateF0 {
public:
    virtual ~EstimateF0() = default;
    virtual bool operator()(F0Contour *output, const Waveform *input) = 0;
};

} }

#endif //UZUME_DSP_ESTIMATE_F0_HPP
