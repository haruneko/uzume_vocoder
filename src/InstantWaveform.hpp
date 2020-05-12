// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_DSP_INSTANT_WAVEFORM_HPP
#define UZUME_DSP_INSTANT_WAVEFORM_HPP

namespace uzume { namespace dsp {

class InstantWaveform {
public:
    const double *wave;
    const double f0;
    const unsigned int length;
    const unsigned int samplingFrequency;
};

} }
#endif //UZUME_DSP_ROOT_INSTANTWAVEFORM_HPP
