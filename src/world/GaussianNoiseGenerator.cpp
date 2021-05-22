// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include "GaussianNoiseGenerator.hpp"

using namespace uzume::dsp::world;

GaussianNoiseGenerator::GaussianNoiseGenerator() {
    reset();
}

void GaussianNoiseGenerator::reset() {
    x = 123456789;
    y = 362436069;
    z = 521288629;
    w = 88675123;
}

double GaussianNoiseGenerator::next() {
    uint32_t t;
    t = x ^ (x << 11);
    x = y;
    y = z;
    z = w;
    w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));

    uint32_t tmp = w >> 4;
    for (int i = 0; i < 11; ++i) {
        t = x ^ (x << 11);
        x = y;
        y = z;
        z = w;
        w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
        tmp += w >> 4;
    }
    return tmp / 268435456.0 - 6.0;
}
