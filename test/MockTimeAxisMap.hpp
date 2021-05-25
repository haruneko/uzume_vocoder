// Copyright 2021 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_VOCODER_MOCK_TIME_AXIS_MAP_HPP
#define UZUME_VOCODER_MOCK_TIME_AXIS_MAP_HPP

#include <gmock/gmock.h>

#include "TimeAxisMap.hpp"

namespace uzume { namespace vocoder {
class MockTimeAxisMap : public TimeAxisMap {
public:
    MOCK_METHOD(double, at, (double), (const, override));
    MOCK_METHOD(double, msLength, (), (const, override));
};
} }

#endif // UZUME_VOCODER_MOCK_TIME_AXIS_MAP_HPP
