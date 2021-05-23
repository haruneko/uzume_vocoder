// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include <cmath>

#include "constant.hpp"
#include "SynthesizeImpulseResponseWithWORLD.hpp"

using namespace uzume::vocoder;
using namespace uzume::vocoder::world;

// WORLD (https://github.com/mmorise/World/) by mmorise is the original of those codes inside namespace.
namespace {

    void fftshift(const double *x, int x_length, double *y) {
        for (int i = 0; i < x_length / 2; ++i) {
            y[i] = x[i + x_length / 2];
            y[i + x_length / 2] = x[i];
        }
    }

    void GetNoiseSpectrum(int noise_size, int fft_size, const ForwardRealFFT *forward_real_fft,
                          GaussianNoiseGenerator *randn) {
        double average = 0.0;
        for (int i = 0; i < noise_size; ++i) {
            forward_real_fft->waveform[i] = randn->next();
            average += forward_real_fft->waveform[i];
        }

        average /= noise_size;
        for (int i = 0; i < noise_size; ++i)
            forward_real_fft->waveform[i] -= average;
        for (int i = noise_size; i < fft_size; ++i)
            forward_real_fft->waveform[i] = 0.0;
        fft_execute(forward_real_fft->forward_fft);
    }

//-----------------------------------------------------------------------------
// GetAperiodicResponse() calculates an aperiodic response.
//-----------------------------------------------------------------------------
    void GetAperiodicResponse(int noise_size, int fft_size,
                              const double *spectrum, const double *aperiodic_ratio, double current_vuv,
                              const ForwardRealFFT *forward_real_fft,
                              const InverseRealFFT *inverse_real_fft,
                              const MinimumPhaseAnalysis *minimum_phase, double *aperiodic_response,
                              GaussianNoiseGenerator *randn) {
        GetNoiseSpectrum(noise_size, fft_size, forward_real_fft, randn);

        if (current_vuv != 0.0)
            for (int i = 0; i <= minimum_phase->fft_size / 2; ++i)
                minimum_phase->log_spectrum[i] = log(spectrum[i] * aperiodic_ratio[i]) / 2.0;
        else
            for (int i = 0; i <= minimum_phase->fft_size / 2; ++i)
                minimum_phase->log_spectrum[i] = log(spectrum[i]) / 2.0;

        GetMinimumPhaseSpectrum(minimum_phase);

        for (int i = 0; i <= fft_size / 2; ++i) {
            inverse_real_fft->spectrum[i][0] =
                    minimum_phase->minimum_phase_spectrum[i][0] * forward_real_fft->spectrum[i][0] -
                    minimum_phase->minimum_phase_spectrum[i][1] * forward_real_fft->spectrum[i][1];
            inverse_real_fft->spectrum[i][1] =
                    minimum_phase->minimum_phase_spectrum[i][0] * forward_real_fft->spectrum[i][1] +
                    minimum_phase->minimum_phase_spectrum[i][1] * forward_real_fft->spectrum[i][0];
        }
        fft_execute(inverse_real_fft->inverse_fft);
        fftshift(inverse_real_fft->waveform, fft_size, aperiodic_response);
    }

//-----------------------------------------------------------------------------
// RemoveDCComponent()
//-----------------------------------------------------------------------------
    void RemoveDCComponent(const double *periodic_response, int fft_size, const double *dc_remover,
                           double *new_periodic_response) {
        double dc_component = 0.0;
        for (int i = fft_size / 2; i < fft_size; ++i)
            dc_component += periodic_response[i];
        for (int i = 0; i < fft_size / 2; ++i)
            new_periodic_response[i] = -dc_component * dc_remover[i];
        for (int i = fft_size / 2; i < fft_size; ++i)
            new_periodic_response[i] -= dc_component * dc_remover[i];
    }

//-----------------------------------------------------------------------------
// GetSpectrumWithFractionalTimeShift() calculates a periodic spectrum with
// the fractional time shift under 1/fs.
//-----------------------------------------------------------------------------
    void GetSpectrumWithFractionalTimeShift(int fft_size, double coefficient, const InverseRealFFT *inverse_real_fft) {
        double re, im, re2, im2;
        for (int i = 0; i <= fft_size / 2; ++i) {
            re = inverse_real_fft->spectrum[i][0];
            im = inverse_real_fft->spectrum[i][1];
            re2 = cos(coefficient * i);
            im2 = sqrt(1.0 - re2 * re2); // sin(pshift)

            inverse_real_fft->spectrum[i][0] = re * re2 + im * im2;
            inverse_real_fft->spectrum[i][1] = im * re2 - re * im2;
        }
    }

//-----------------------------------------------------------------------------
// GetPeriodicResponse() calculates a periodic response.
//-----------------------------------------------------------------------------
    void GetPeriodicResponse(int fft_size, const double *spectrum,
                             const double *aperiodic_ratio, double current_vuv,
                             const InverseRealFFT *inverse_real_fft,
                             const MinimumPhaseAnalysis *minimum_phase, const double *dc_remover,
                             double fractional_time_shift, int fs, double *periodic_response) {
        if (current_vuv <= 0.5 || aperiodic_ratio[0] > 0.999) {
            for (int i = 0; i < fft_size; ++i)
                periodic_response[i] = 0.0;
            return;
        }

        for (int i = 0; i <= minimum_phase->fft_size / 2; ++i)
            minimum_phase->log_spectrum[i] = log(spectrum[i] * (1.0 - aperiodic_ratio[i]) + SafeGuardMinimum) / 2.0;
        GetMinimumPhaseSpectrum(minimum_phase);

        for (int i = 0; i <= fft_size / 2; ++i) {
            inverse_real_fft->spectrum[i][0] = minimum_phase->minimum_phase_spectrum[i][0];
            inverse_real_fft->spectrum[i][1] = minimum_phase->minimum_phase_spectrum[i][1];
        }

        // apply fractional time delay of fractional_time_shift seconds
        // using linear phase shift
        double coefficient = 2.0 * Pi * fractional_time_shift * fs / fft_size;
        GetSpectrumWithFractionalTimeShift(fft_size, coefficient, inverse_real_fft);

        fft_execute(inverse_real_fft->inverse_fft);
        fftshift(inverse_real_fft->waveform, fft_size, periodic_response);
        RemoveDCComponent(periodic_response, fft_size, dc_remover, periodic_response);
    }

} // namespace

SynthesizeImpulseResponseWithWORLD::SynthesizeImpulseResponseWithWORLD(unsigned int fftSize, unsigned int samplingFrequency)
        : fftSize(fftSize), samplingFrequency(samplingFrequency),
          forwardRealFFT({0}), inverseRealFFT({0}), minimumPhaseAnalysis({0}),
          dcRemover(nullptr), periodicResponse(nullptr), aperiodicResponse(nullptr), randn() {
    if (fftSize != 0 && samplingFrequency != 0) {
        InitializeForwardRealFFT(fftSize, &forwardRealFFT);
        InitializeInverseRealFFT(fftSize, &inverseRealFFT);
        InitializeMinimumPhaseAnalysis(fftSize, &minimumPhaseAnalysis);
        dcRemover = new double[fftSize];
        double dcComponent = 0.0;
        for (unsigned int i = 0; i < fftSize / 2; ++i) {
            dcRemover[i] = 0.5 - 0.5 * cos(2.0 * Pi * (i + 1.0) / (1.0 + fftSize));
            dcRemover[fftSize - i - 1] = dcRemover[i];
            dcComponent += dcRemover[i] * 2.0;
        }
        for (unsigned int i = 0; i < fftSize / 2; ++i) {
            dcRemover[i] /= dcComponent;
            dcRemover[fftSize - i - 1] = dcRemover[i];
        }
        periodicResponse = new double[fftSize];
        aperiodicResponse = new double[fftSize];
    }
}

ImpulseResponseParameters::ImpulseResponseParameters(unsigned int fftSize)
        : spectrum(nullptr), secondsFractionalTimeShift(0.0), VUV(0.0), noiseSize(0) {
    if (fftSize != 0) {
        spectrum = new Spectrum(fftSize);
    }
}

ImpulseResponseParameters::~ImpulseResponseParameters() {
    delete spectrum;
}

SynthesizeImpulseResponseWithWORLD::~SynthesizeImpulseResponseWithWORLD() {
    DestroyForwardRealFFT(&forwardRealFFT);
    DestroyInverseRealFFT(&inverseRealFFT);
    DestroyMinimumPhaseAnalysis(&minimumPhaseAnalysis);
    delete[] dcRemover;
    delete[] periodicResponse;
    delete[] aperiodicResponse;
}

bool SynthesizeImpulseResponseWithWORLD::operator()(ImpulseResponse *output, const ImpulseResponseParameters *frame) {
    if (output->length() < frame->spectrum->fftSize) {
        return false;
    }
    // Synthesis of the periodic response
    GetPeriodicResponse(fftSize, frame->spectrum->periodicSpectrum, frame->spectrum->aperiodicSpectrum,
                        frame->VUV, &inverseRealFFT, &minimumPhaseAnalysis, dcRemover,
                        frame->secondsFractionalTimeShift, samplingFrequency, periodicResponse);

    // Synthesis of the aperiodic response
    GetAperiodicResponse(frame->noiseSize, fftSize, frame->spectrum->periodicSpectrum,
                         frame->spectrum->aperiodicSpectrum, frame->VUV, &forwardRealFFT,
                         &inverseRealFFT, &minimumPhaseAnalysis, aperiodicResponse, &randn);

    double sqrt_noise_size = sqrt(static_cast<double>(frame->noiseSize));
    double *response = output->data();
    for (unsigned int i = 0; i < fftSize; ++i) {
        response[i] = (periodicResponse[i] * sqrt_noise_size + aperiodicResponse[i]) / fftSize;
    }
    return true;
}
