// Copyright 2021 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_VOCODER_TIME_AXIS_MAP_HPP
#define UZUME_VOCODER_TIME_AXIS_MAP_HPP

namespace uzume { namespace vocoder {

/**
 * TimeAxisMap represents a time axis map.
 * This class will be used to stretch Spectrogram, Contour or other time based values.
 */
class TimeAxisMap {
public:
    TimeAxisMap() = default;
    virtual ~TimeAxisMap() = default;
    /**
     * When m is a TimeAxisMap, m will be used like `g(t) = f(m(t))`.
     * @param ms represents in `g` function.
     * @return the corresponding time in `f` function.
     */
    virtual double at(double ms) const = 0;

    /**
     * @return an entire length of this map in milli seconds.
     */
    virtual double msLength() const = 0;
};

} }
#endif //UZUME_VOCODER_TIME_AXIS_MAP_HPP
