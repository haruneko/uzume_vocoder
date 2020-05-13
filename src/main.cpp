// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.#include <stdio.h>
#include <algorithm>
#include <cmath>

#include "constant.hpp"
#include "AnalyzeAperiodicityWithD4C.hpp"
#include "AnalyzePeriodicityWithCheapTrick.hpp"
#include "EstimateF0WithDIO.hpp"
#include "NaiveSpectrogram.hpp"
#include "SynthesizeImpulseResponseWithWORLD.hpp"
#include "SynthesizePhraseWithWORLD.hpp"

using namespace uzume::dsp;

// implement wave read, write function.
int GetAudioLength(const char*);
void wavread(const char*, int*, int*, double*);
void wavwrite(const double*, int, int, int, const char*);

// This is a sample to use dsp directory.
int main()
{

    const char *inputPath = "/path/to/input.wav";
    const char *outputPath = "/path/to/output.wav";
    Waveform waveform;
    waveform.length = GetAudioLength(inputPath);
    waveform.data = new double[waveform.length];
    int waveLength, samplingFrequency, samplingBits;
    wavread(inputPath, &samplingFrequency, &samplingBits, waveform.data);
    waveform.samplingFrequency = samplingFrequency;

    int fftSize = pow(2.0, 1.0 + (int)(log(3.0 * waveform.samplingFrequency / 71.0 + 1) / Log2));

    /* Analyze and create WORLD spectrogram here. */
    EstimateF0WithDIO dio(2.0);
    F0Contour f0;
    f0.length = dio.getF0LengthFor(waveform.samplingFrequency, waveform.length);
    f0.msFramePeriod = 2.0;
    f0.data = new double[f0.length];
    dio(&f0, &waveform);

    NaiveSpectrogram spectrogram((unsigned int)f0.length, (unsigned int)fftSize, f0.msFramePeriod);

    for(int i = 0; i < f0.length; i++) {
        spectrogram.f0Contour[i] = f0.data[i];
    }

    AnalyzeAperiodicityWithD4C d4c(fftSize, samplingFrequency);
    AnalyzePeriodicityWithCheapTrick cheapTrick(fftSize, samplingFrequency);
    InstantWaveform iw;
    Spectrum sp(fftSize);
    iw.wave = new double[fftSize * 2] + fftSize;
    iw.samplingFrequency = waveform.samplingFrequency;
    iw.length = fftSize * 2;
    for(int i = 0; i < f0.length; i++) {
        int index = i * f0.msFramePeriod / 1000.0 * waveform.samplingFrequency;
        iw.f0 = f0.data[i];
        for(int wi = -fftSize; wi < fftSize; wi++) {
            iw.wave[wi] = (index + wi < 0) ? 0.0 :
                                    (index + wi >= waveform.length) ? 0.0 :
                                    waveform.data[index + wi];
        }
        cheapTrick(&sp, &iw);
        d4c(&sp, &iw);
        for(int j = 0; j <= sp.fftSize / 2; j++) {
            spectrogram.periodicSpecgram[i][j] = sp.periodicSpectrum[j];
            spectrogram.aperiodicSpecgram[i][j] = sp.aperiodicSpectrum[j];
        }
    }
    delete[] (iw.wave - fftSize);

    for(int i = 0; i < waveform.length; i++) {
        waveform.data[i] = 0.0;
    }
    SynthesizeImpulseResponseWithWORLD irs(spectrogram.fftSize(), samplingFrequency);
    SynthesizePhraseWithWORLD synthesize(&irs);

    PhraseSignal s(waveform.data, /* indexMin = */ 0, /* indexMax = */ waveform.length, samplingFrequency);
    PhraseParameters p(&spectrogram, /* startPhase = */ 0.0, /* startFractionalTimeShift = */ 0.0);

    synthesize(&s, &p);

    /* Save wave as wav file here. */
    wavwrite(waveform.data, waveform.length, waveform.samplingFrequency, 16, outputPath);
}
