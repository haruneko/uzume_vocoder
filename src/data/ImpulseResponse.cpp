// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include "ImpulseResponse.hpp"

using namespace uzume::vocoder;

ImpulseResponse::ImpulseResponse(unsigned int length) : raw(nullptr), _length(length) {
    if (_length != 0) {
        raw = new double[_length];
    }
}

ImpulseResponse::~ImpulseResponse() {
    delete[] raw;
}
