// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_DSP_SYNTHESIZE_IMPULSE_RESPONSE_WITH_WORLD_HPP
#define UZUME_DSP_SYNTHESIZE_IMPULSE_RESPONSE_WITH_WORLD_HPP
#include "fft.hpp"
#include "GaussianNoiseGenerator.hpp"
#include "SynthesizeImpulseResponse.hpp"

namespace uzume { namespace dsp

/**
* SynthesizeImpulseResponseWithWORLD is WORLD's implementation of the SynthesizeImpulseResponse class.
*/
class SynthesizeImpulseResponseWithWORLD final : public SynthesizeImpulseResponse {
public:
    SynthesizeImpulseResponseWithWORLD(unsigned int fftSize, unsigned int samplingFrequency);

    ~SynthesizeImpulseResponseWithWORLD() override;

    /**
     * SynthesizeImpulseResponseWithWORLD#() synthesizes a single impulse response from WORLD's spectrum and aperiodicity.
     * @param output is an output signal of synthesis. output should have a buffer enough long for FFT size of synthesis.
     * @param frame is a single frame that contains WORLD's spectral envelope and aperiodic ratio.
     * @return true if synthesis succeeds.
     *         false if synthesis fails.
     */
    bool operator()(ImpulseResponse *output, const ImpulseResponseParameters *input) override;

    const unsigned int fftSize;
    const unsigned int samplingFrequency;
private:
    // Those below are used as buffer...
    ForwardRealFFT forwardRealFFT;
    InverseRealFFT inverseRealFFT;
    MinimumPhaseAnalysis minimumPhaseAnalysis;
    double *dcRemover;

    double *periodicResponse;
    double *aperiodicResponse;

    GaussianNoiseGenerator randn;
};

} }
#endif //UZUME_DSP_SYNTHESIZE_IMPULSE_RESPONSE_WITH_WORLD_HPP
