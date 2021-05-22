// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_DSP_CONSTANT_HPP
#define UZUME_DSP_CONSTANT_HPP

namespace uzume { namespace dsp { namespace world {
    extern const double CutOff;
    extern const double FloorF0StoneMask;
    extern const double Pi;
    extern const double SafeGuardMinimum;
    extern const double Eps;
    extern const double FloorF0;
    extern const double CeilF0;
    extern const double DefaultF0;
    extern const double Log2;
    extern const double MaximumValue;
    extern const double DefaultQ1;
    enum WindowFunctionType {
        Hanning = 1,
        Blackman,
    };
    extern const double FrequencyInterval;
    extern const double UpperLimit;
    extern const double Threshold;
    extern const double FloorF0D4C;
    extern const double LowestF0D4CLoveTrain;
} } }

#endif //UZUME_DSP_CONSTANT_HPP
