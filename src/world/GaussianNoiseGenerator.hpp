// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef __GAUSSIANNOISEGENERATOR_HPP__
#define __GAUSSIANNOISEGENERATOR_HPP__

#include <stdint.h>

namespace uzume { namespace dsp { namespace world {

/**
 * GaussianNoiseGenerator is a simple C++ implementation of `randn` function in WORLD (https://github.com/mmorise/World/).
 */
class GaussianNoiseGenerator final {
public:
    GaussianNoiseGenerator();

    /**
     * next returns randomized number. This method changes generator's state, so this function is a thread-unsafe one.
     */
    double next();

    void reset();

private:
    uint32_t x;
    uint32_t y;
    uint32_t z;
    uint32_t w;
};

} } }

#endif
