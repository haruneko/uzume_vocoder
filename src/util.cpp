// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include <cmath>

#include "constant.hpp"
#include "util.hpp"

using namespace uzume::dsp;

namespace {
    void diff(const double *x, int x_length, double *y) {
        for (int i = 0; i < x_length - 1; ++i) y[i] = x[i + 1] - x[i];
    }

    void interp1Q(double x, double shift, const double *y, int x_length,
                  const double *xi, int xi_length, double *yi) {
        double *xi_fraction = new double[xi_length];
        double *delta_y = new double[x_length];
        int *xi_base = new int[xi_length];

        double delta_x = shift;
        for (int i = 0; i < xi_length; ++i) {
            xi_base[i] = static_cast<int>((xi[i] - x) / delta_x);
            xi_fraction[i] = (xi[i] - x) / delta_x - xi_base[i];
        }
        diff(y, x_length, delta_y);
        delta_y[x_length - 1] = 0.0;

        for (int i = 0; i < xi_length; ++i)
            yi[i] = y[xi_base[i]] + delta_y[xi_base[i]] * xi_fraction[i];

        // Bug was fixed at 2013/07/14 by M. Morise
        delete[] xi_fraction;
        delete[] xi_base;
        delete[] delta_y;
    }
    void SetParametersForLinearSmoothing(int boundary, int fft_size, int fs,
                                         double width, const double *power_spectrum, double *mirroring_spectrum,
                                         double *mirroring_segment, double *frequency_axis) {
        for (int i = 0; i < boundary; ++i)
            mirroring_spectrum[i] = power_spectrum[boundary - i];
        for (int i = boundary; i < fft_size / 2 + boundary; ++i)
            mirroring_spectrum[i] = power_spectrum[i - boundary];
        for (int i = fft_size / 2 + boundary; i <= fft_size / 2 + boundary * 2; ++i)
            mirroring_spectrum[i] = power_spectrum[fft_size / 2 - (i - (fft_size / 2 + boundary))];

        mirroring_segment[0] = mirroring_spectrum[0] * fs / fft_size;
        for (int i = 1; i < fft_size / 2 + boundary * 2 + 1; ++i)
            mirroring_segment[i] = mirroring_spectrum[i] * fs / fft_size + mirroring_segment[i - 1];

        for (int i = 0; i <= fft_size / 2; ++i)
            frequency_axis[i] = static_cast<double>(i) / fft_size * fs - width / 2.0;
    }

}

int uzume::dsp::matlab_round(double x) {
    return x > 0 ? static_cast<int>(x + 0.5) : static_cast<int>(x - 0.5);
}


void uzume::dsp::DCCorrection(const double *input, double f0, int fs, int fft_size, double *output) {
    int upper_limit = 2 + static_cast<int>(f0 * fft_size / fs);
    double *low_frequency_replica = new double[upper_limit];
    double *low_frequency_axis = new double[upper_limit];

    for (int i = 0; i < upper_limit; ++i)
        low_frequency_axis[i] = static_cast<double>(i) * fs / fft_size;

    int upper_limit_replica = upper_limit - 1;
    interp1Q(f0 - low_frequency_axis[0],
             -static_cast<double>(fs) / fft_size, input, upper_limit + 1,
             low_frequency_axis, upper_limit_replica, low_frequency_replica);

    for (int i = 0; i < upper_limit_replica; ++i)
        output[i] = input[i] + low_frequency_replica[i];

    delete[] low_frequency_replica;
    delete[] low_frequency_axis;
}

void uzume::dsp::LinearSmoothing(const double *input, double width, int fs, int fft_size, double *output) {
    int boundary = static_cast<int>(width * fft_size / fs) + 1;

    // These parameters are set by the other function.
    double *mirroring_spectrum = new double[fft_size / 2 + boundary * 2 + 1];
    double *mirroring_segment = new double[fft_size / 2 + boundary * 2 + 1];
    double *frequency_axis = new double[fft_size / 2 + 1];
    SetParametersForLinearSmoothing(boundary, fft_size, fs, width, input, mirroring_spectrum, mirroring_segment,
                                    frequency_axis);

    double *low_levels = new double[fft_size / 2 + 1];
    double *high_levels = new double[fft_size / 2 + 1];
    double origin_of_mirroring_axis = -(boundary - 0.5) * fs / fft_size;
    double discrete_frequency_interval = static_cast<double>(fs) / fft_size;

    interp1Q(origin_of_mirroring_axis, discrete_frequency_interval,
             mirroring_segment, fft_size / 2 + boundary * 2 + 1, frequency_axis,
             fft_size / 2 + 1, low_levels);

    for (int i = 0; i <= fft_size / 2; ++i) frequency_axis[i] += width;

    interp1Q(origin_of_mirroring_axis, discrete_frequency_interval,
             mirroring_segment, fft_size / 2 + boundary * 2 + 1, frequency_axis,
             fft_size / 2 + 1, high_levels);

    for (int i = 0; i <= fft_size / 2; ++i)
        output[i] = (high_levels[i] - low_levels[i]) / width;

    delete[] mirroring_spectrum;
    delete[] mirroring_segment;
    delete[] frequency_axis;
    delete[] low_levels;
    delete[] high_levels;
}