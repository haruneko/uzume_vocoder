// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_VOCODER_ESTIMATE_F0_HPP
#define UZUME_VOCODER_ESTIMATE_F0_HPP
#include "data/Contour.hpp"
#include "data/Waveform.hpp"

namespace uzume { namespace vocoder {

/**
 * EstimateF0 is an interface to estimate f0.
 */
class EstimateF0 {
public:
    virtual ~EstimateF0() = default;
    /**
     * () estimates f0 contour from input and put it into output.
     * @return true if estimation is successful, otherwise false.
     */
    virtual bool operator()(Contour *output, const Waveform *input) = 0;
};

} }

#endif //UZUME_VOCODER_ESTIMATE_F0_HPP
