// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef __DSP_SYNTHESIZE_IMPULSE_RESPONSE_HPP__
#define __DSP_SYNTHESIZE_IMPULSE_RESPONSE_HPP__

#include "fft.hpp"
#include "GaussianNoiseGenerator.hpp"
#include "ImpulseResponse.hpp"
#include "Spectrum.hpp"

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

#endif
