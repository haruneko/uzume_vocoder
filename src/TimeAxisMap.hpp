// Copyright 2021 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_VOCODER_TIME_AXIS_MAP_HPP
#define UZUME_VOCODER_TIME_AXIS_MAP_HPP

namespace uzume { namespace vocoder {

class TimeAxisMap {
public:
    TimeAxisMap() = default;
    virtual ~TimeAxisMap() = default;
    virtual double at(double ms) const = 0;
    virtual double msLength() const = 0;
};

} }
#endif //UZUME_VOCODER_TIME_AXIS_MAP_HPP
