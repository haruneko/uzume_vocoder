// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_DSP_UTIL_HPP
#define UZUME_DSP_UTIL_HPP

#include "fft.hpp"
#include "GaussianNoiseGenerator.hpp"

namespace uzume { namespace dsp {
template<class T>
inline T myMin(T x, T y) {
    return (x < y) ? x : y;
}

template<class T>
inline T myMax(T x, T y) {
    return (x > y) ? x : y;
}

int matlab_round(double x);

void DCCorrection(const double *input, double f0, int fs, int fft_size, double *output);

void LinearSmoothing(const double *input, double width, int fs, int fft_size, double *output);

} }

#endif //UZUME_DSP_UTIL_HPP
