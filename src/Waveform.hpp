// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_DSP_WAVEFORM_HPP
#define UZUME_DSP_WAVEFORM_HPP

namespace uzume { namespace dsp {

class Waveform {
public:
    double *data;
    unsigned int length;
    unsigned int samplingFrequency;
};

} }

#endif //UZUME_DSP_ROOT_WAVEFORM_HPP
