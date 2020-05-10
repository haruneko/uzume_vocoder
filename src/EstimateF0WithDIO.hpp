// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_DSP_ESTIMATE_F0_WITH_DIO_HPP
#define UZUME_DSP_ESTIMATE_F0_WITH_DIO_HPP
#include "EstimateF0.hpp"

namespace uzume { namespace dsp {

class EstimateF0WithDio : public EstimateF0 {
public:
    EstimateF0WithDio() = delete;
    EstimateF0WithDio(double msFramePeriod);

    bool operator()(F0Contour *output, const Waveform *input) override;
    int getF0LengthFor(unsigned int samplingFrequency, unsigned int waveLength) const;
private:
    double msFramePeriod;
};

} }

#endif //UZUME_DSP_ESTIMATE_F0_WITH_DIO_HPP