// Modified by Hal@shurabaP.
// The original of this source code is WORLD(https://github.com/mmorise/World/) by mmorise.
#ifndef __DSP_FFT_HPP__
#define __DSP_FFT_HPP__

#include "fft_util.hpp"

fft_plan fft_plan_dft_1d(int n, fft_complex *in, fft_complex *out, int sign, unsigned int flags);

fft_plan fft_plan_dft_c2r_1d(int n, fft_complex *in, double *out, unsigned int flags);

fft_plan fft_plan_dft_r2c_1d(int n, double *in, fft_complex *out, unsigned int flags);

void fft_execute(fft_plan p);

void fft_destroy_plan(fft_plan p);

#endif
