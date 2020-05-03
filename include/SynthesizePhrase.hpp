// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef __DSP_SYNTHESIZE_PHRASE_HPP__
#define __DSP_SYNTHESIZE_PHRASE_HPP__

#include "Spectrogram.hpp"
#include "SynthesizeImpulseResponse.hpp"

/**
 * PhraseSignal is a value class for phrase synthesis output.
 * This class is not responsible for memory allocation.
 * Make sure neither PhraseSignal#indexMin nor PhraseSignal#indexMax is out of array range.
 */
class PhraseSignal final {
public:
    PhraseSignal() = delete;
    PhraseSignal(double *waveStartPoint, int indexMin, int indexMax, unsigned int samplingFrequency);
    double *waveStartPoint;
    int indexMin;
    int indexMax;
    unsigned int samplingFrequency;
};

/**
 * PhraseSignal is a value class for phrase synthesis input.
 * This class is not responsible for memory allocation.
 * PhraseParameters#startPhase and PhraseParameters#startFractionalTimeShift should be set
 * when SynthesizePhrase begins synthesis from the middle of music sentence, otherwise they should be 0.0.
 */
class PhraseParameters final {
public:
    PhraseParameters() = delete;
    PhraseParameters(Spectrogram *spectrogram, double startPhase, double startFractionalTimeShift);

    Spectrogram *spectrogram;
    double startPhase;
    double startFractionalTimeShift;
};

/**
 * SynthesizePhrase synthesize a single phrase.
 */
class SynthesizePhrase {
public:
    SynthesizePhrase() = default;
    SynthesizePhrase(const SynthesizePhrase&) = delete;
    virtual ~SynthesizePhrase() = default;

    /**
     * SynthesizePhrase#() synthesize a single phrase into output from input.
     * @return true if synthesis is successful.
     * @return false otherwise.
     */
    virtual bool operator()(PhraseSignal *output, PhraseParameters *input) = 0;
};

class SynthesizePhraseWithWORLD final : public SynthesizePhrase {
public:
    SynthesizePhraseWithWORLD() = delete;
    explicit SynthesizePhraseWithWORLD(SynthesizeImpulseResponse *synthesize, double f0Floor = 71.0, double f0Default = 500.0);
    ~SynthesizePhraseWithWORLD() override = default;

    bool operator()(PhraseSignal *output, PhraseParameters *input) override;
private:
    SynthesizeImpulseResponse *synthesize;

    double f0Floor;
    double f0Default;
};

#endif
