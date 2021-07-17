// Copyright 2021 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include "LinearTimeAxisMap.hpp"

using namespace uzume::vocoder;

LinearTimeAxisMap::LinearTimeAxisMap(double msLeft, double msRight, double msLength)
    : msLeft(msLeft), msRight(msRight), _msLength(msLength) {
}

double LinearTimeAxisMap::at(double ms) const {
    return msLeft + ms / _msLength * (msRight - msLeft);
}

double LinearTimeAxisMap::msLength() const {
    return _msLength;
}
