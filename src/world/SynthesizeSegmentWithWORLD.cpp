// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include <cmath>
#include <algorithm>

#include "constant.hpp"
#include "SynthesizeSegmentWithWORLD.hpp"

using namespace uzume::dsp::world;

SynthesizeSegmentWithWORLD::SynthesizeSegmentWithWORLD(SynthesizeImpulseResponse *synthesize)
        : synthesize(synthesize) {
}

bool SynthesizeSegmentWithWORLD::operator()(SegmentSignal *output, SegmentParameters *input) {

    ImpulseResponse responseBuffer(input->spectrogram->fftSize());
    ImpulseResponseParameters frameBuffer(input->spectrogram->fftSize());

    double *wave = output->waveStartPoint;
    unsigned int samplingFrequency = output->samplingFrequency;
    auto waveLength = (unsigned int)(input->spectrogram->msLength() / 1000.0 * samplingFrequency);

    unsigned int fftSize = input->spectrogram->fftSize();

    double previousPhase = input->startPhase;
    double currentPhase = input->startPhase;
    double previousFractionalTimeShift = input->startFractionalTimeShift;
    double previousVUV = input->spectrogram->f0At(0.0) < FloorF0? 0.0: 1.0;
    int previousPulseIndex = 0;

    for(unsigned int i = 0; i < waveLength; i++) {
        double ms = (double)i / samplingFrequency * 1000.0;
        double f0Interpolated = input->spectrogram->f0At(ms);
        double currentVUV = f0Interpolated < FloorF0? 0.0 : 1.0;
        double currentF0 = f0Interpolated < FloorF0? DefaultF0 : f0Interpolated;
        currentPhase += 2 * Pi * currentF0 / samplingFrequency;
        currentPhase = fmod(currentPhase, 2 * Pi);

        // Pulse occurs.
        if(fabs(currentPhase - previousPhase) > Pi) {
            double y1 = currentPhase - 2.0 * Pi;
            double y2 = previousPhase;
            double x = -y1 / (y2 - y1);
            double currentFractionalTimeShift = x / samplingFrequency;

            frameBuffer.VUV = previousVUV;
            frameBuffer.noiseSize = i - previousPulseIndex;
            frameBuffer.secondsFractionalTimeShift = previousFractionalTimeShift;
            input->spectrogram->pickUpSpectrumAt(frameBuffer.spectrum, ms);
            // calculate safe aperiodicity ratio.
            for(unsigned int i = 0; i <= fftSize / 2; i++) {
                double safeValue = std::max<double>(0.001, std::min<double>(0.999999999999, frameBuffer.spectrum->aperiodicSpectrum[i]));
                frameBuffer.spectrum->aperiodicSpectrum[i] = pow(safeValue, 2.0);
            }

            (*synthesize)(&responseBuffer, &frameBuffer);
            int waveIndexOffset = previousPulseIndex - fftSize / 2 + 1;
            double *impulseResponse = responseBuffer.data();
            int begin = std::max(std::max(-waveIndexOffset + output->indexMin, 0), -waveIndexOffset);

            int end = std::min(std::min(output->indexMax - waveIndexOffset, (int)fftSize), (int)waveLength - waveIndexOffset);
            for (int responseIndex = begin; responseIndex < end; ++responseIndex) {
                int waveIndex = responseIndex + waveIndexOffset;
                wave[waveIndex] += impulseResponse[responseIndex];
            }

            previousPulseIndex = i;
            previousFractionalTimeShift = currentFractionalTimeShift;
            previousVUV = currentVUV;
        }
        previousPhase = currentPhase;
    }
    return true;
}
