// Copyright 2021 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef __UZUME_VOCODER_LINEAR_TIME_AXIS_MAP__
#define __UZUME_VOCODER_LINEAR_TIME_AXIS_MAP__

#include "../TimeAxisMap.hpp"

namespace uzume {
namespace vocoder {

/**
 * LinearTimeAxisMap is an implementation of TimeAxisMap that simply cut and stretch time axis.
 */
class LinearTimeAxisMap final : public TimeAxisMap {
public:
    LinearTimeAxisMap(double msLeft, double msRight, double msLength);
    double at(double ms) const override;
    double msLength() const override;
private:
    double msLeft;
    double msRight;
    double _msLength;
};

}
}

#endif // __UZUME_VOCODER_LINEAR_TIME_AXIS_MAP__
