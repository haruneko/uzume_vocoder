// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef __DSP_SYNTHESIZE_IMPULSE_RESPONSE_HPP__
#define __DSP_SYNTHESIZE_IMPULSE_RESPONSE_HPP__

#include "data/ImpulseResponse.hpp"
#include "data/Spectrum.hpp"

namespace uzume { namespace dsp {

class ImpulseResponseParameters final {
public:
    explicit ImpulseResponseParameters(unsigned int fftSize);

    ~ImpulseResponseParameters();

    Spectrum *spectrum;
    // Pulse is located on the position; (index / sampling frequency)[second] + secondsFractionalTimeShift[second].
    double secondsFractionalTimeShift;
    double VUV;
    int noiseSize;
};

/**
 * SynthesizeImpulseResponse is a class for the function object to synthesize a single impulse response.
 */
class SynthesizeImpulseResponse {
public:
    SynthesizeImpulseResponse &operator=(const SynthesizeImpulseResponse&) = delete;
    virtual ~SynthesizeImpulseResponse() = default;
    /**
     * SynthesizeImpulseResponse#() synthesizes a single impulse response.
     * @param output is an output signal of synthesis. output should have a buffer enough long for FFT size of synthesis.
     * @param frame is a single frame to synthesize.
     * @return true if synthesis succeeds.
     *         false if synthesis fails.
     */
    virtual bool operator()(ImpulseResponse *output, const ImpulseResponseParameters *input) = 0;
};


} }

#endif
