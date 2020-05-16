// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include <algorithm>
#include <cmath>
#include "Contour.hpp"

using namespace uzume::dsp;

Contour::Contour(double msLength, double msFramePeriod)
        : length((int)(msLength / msFramePeriod) + 1), data(nullptr), msFramePeriod(msFramePeriod) {

    if(length > 0) {
        data = new double[length];
    }
}

Contour::~Contour() {
    delete[] data;
}

double Contour::at(double ms) const {
    int il = (int)floor(ms / msFramePeriod);
    int ir = il + 1;
    double interpolation = ms / msFramePeriod - (double)il;
    double vl = data[std::max(0, std::min(length - 1, il))];
    double vr = data[std::max(0, std::min(length - 1, ir))];

    return vl * (1 - interpolation) + vr * interpolation;
}

double Contour::msLength() const {
    return msFramePeriod * length;
}
