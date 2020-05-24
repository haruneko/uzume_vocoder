// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.#include <stdio.h>
#include <algorithm>

#include "SynthesizeImpulseResponseWithWORLD.hpp"
#include "SynthesizePhraseWithWORLD.hpp"
#include "Waveform.hpp"
#include "WaveformSpectrogram.hpp"

using namespace uzume::dsp;

// This is a sample to use dsp directory.
int main() {
    const char *inputPath = "/path/to/input.wav";
    const char *outputPath = "/path/to/output.wav";
    Waveform *input = Waveform::read(inputPath);
    Waveform *output = new Waveform(input->length, input->samplingFrequency);
    WaveformSpectrogram spectrogram(input);

    SynthesizeImpulseResponseWithWORLD irs(spectrogram.fftSize(), input->samplingFrequency);
    SynthesizePhraseWithWORLD synthesize(&irs);

    for(unsigned int i = 0; i < output->length; i++) {
        output->data[i] = 0.0;
    }

    PhraseSignal s(output->data, /* indexMin = */ 0, /* indexMax = */ output->length, output->samplingFrequency);
    PhraseParameters p(&spectrogram, /* startPhase = */ 0.0, /* startFractionalTimeShift = */ 0.0);

    synthesize(&s, &p);

    output->save(outputPath);
}
