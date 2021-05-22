// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include "AnalyzePeriodicityWithCheapTrick.hpp"

#include <cmath>

#include "constant.hpp"
#include "util.hpp"

using namespace uzume::dsp::world;

namespace {

    //-----------------------------------------------------------------------------
    // SmoothingWithRecovery() carries out the spectral smoothing and spectral
    // recovery on the Cepstrum domain.
    //-----------------------------------------------------------------------------
    void SmoothingWithRecovery(double f0, int fs, int fft_size, double q1,
                               const ForwardRealFFT *forward_real_fft,
                               const InverseRealFFT *inverse_real_fft, double *spectral_envelope) {
        double *smoothing_lifter = new double[fft_size];
        double *compensation_lifter = new double[fft_size];

        smoothing_lifter[0] = 1.0;
        compensation_lifter[0] = (1.0 - 2.0 * q1) + 2.0 * q1;
        double quefrency;
        for (int i = 1; i <= forward_real_fft->fft_size / 2; ++i) {
            quefrency = static_cast<double>(i) / fs;
            smoothing_lifter[i] = sin(Pi * f0 * quefrency) / (Pi * f0 * quefrency);
            compensation_lifter[i] = (1.0 - 2.0 * q1) + 2.0 * q1 * cos(2.0 * Pi * quefrency * f0);
        }

        for (int i = 0; i <= fft_size / 2; ++i)
            forward_real_fft->waveform[i] = log(forward_real_fft->waveform[i]);
        for (int i = 1; i < fft_size / 2; ++i)
            forward_real_fft->waveform[fft_size - i] = forward_real_fft->waveform[i];
        fft_execute(forward_real_fft->forward_fft);

        for (int i = 0; i <= fft_size / 2; ++i) {
            inverse_real_fft->spectrum[i][0] = forward_real_fft->spectrum[i][0] *
                                               smoothing_lifter[i] * compensation_lifter[i] / fft_size;
            inverse_real_fft->spectrum[i][1] = 0.0;
        }
        fft_execute(inverse_real_fft->inverse_fft);

        for (int i = 0; i <= fft_size / 2; ++i)
            spectral_envelope[i] = exp(inverse_real_fft->waveform[i]);

        delete[] smoothing_lifter;
        delete[] compensation_lifter;
    }

    //-----------------------------------------------------------------------------
    // GetPowerSpectrum() calculates the power_spectrum with DC correction.
    // DC stands for Direct Current. In this case, the component from 0 to F0 Hz
    // is corrected.
    //-----------------------------------------------------------------------------
    void GetPowerSpectrum(int fs, double f0, int fft_size,
                          const ForwardRealFFT *forward_real_fft) {
        int half_window_length = matlab_round(1.5 * fs / f0);

        // FFT
        for (int i = half_window_length * 2 + 1; i < fft_size; ++i)
            forward_real_fft->waveform[i] = 0.0;
        fft_execute(forward_real_fft->forward_fft);

        // Calculation of the power spectrum.
        double *power_spectrum = forward_real_fft->waveform;
        for (int i = 0; i <= fft_size / 2; ++i)
            power_spectrum[i] =
                    forward_real_fft->spectrum[i][0] * forward_real_fft->spectrum[i][0] +
                    forward_real_fft->spectrum[i][1] * forward_real_fft->spectrum[i][1];

        // DC correction
        DCCorrection(power_spectrum, f0, fs, fft_size, power_spectrum);
    }

    //-----------------------------------------------------------------------------
    // AddInfinitesimalNoise()
    //-----------------------------------------------------------------------------
    void AddInfinitesimalNoise(const double *input_spectrum, int fft_size,
                               double *output_spectrum, GaussianNoiseGenerator *randn) {
        for (int i = 0; i <= fft_size / 2; ++i)
            output_spectrum[i] = input_spectrum[i] + fabs(randn->next()) * Eps;
    }

    //-----------------------------------------------------------------------------
    // SetParametersForGetWindowedWaveform()
    //-----------------------------------------------------------------------------
    void SetParametersForGetWindowedWaveform(int half_window_length, int fs, double current_f0, int *base_index,
                                             double *window) {
        for (int i = -half_window_length; i <= half_window_length; ++i)
            base_index[i + half_window_length] = i;

        // Designing of the window function
        double average = 0.0;
        double position;
        for (int i = 0; i <= half_window_length * 2; ++i) {
            position = base_index[i] / 1.5 / fs;
            window[i] = 0.5 * cos(Pi * position * current_f0) + 0.5;
            average += window[i] * window[i];
        }
        average = sqrt(average);
        for (int i = 0; i <= half_window_length * 2; ++i) window[i] /= average;
    }

    //-----------------------------------------------------------------------------
    // GetWindowedWaveform() windows the waveform by F0-adaptive window
    //-----------------------------------------------------------------------------
    void GetWindowedWaveform(const double *x, int fs, double current_f0, const ForwardRealFFT *forward_real_fft, GaussianNoiseGenerator *randn) {
        int half_window_length = matlab_round(1.5 * fs / current_f0);

        int *base_index = new int[half_window_length * 2 + 1];
        double *window = new double[half_window_length * 2 + 1];

        SetParametersForGetWindowedWaveform(half_window_length, fs, current_f0, base_index, window);

        // F0-adaptive windowing
        double *waveform = forward_real_fft->waveform;
        for (int i = 0; i <= half_window_length * 2; ++i)
            waveform[i] = x[base_index[i]] * window[i] + randn->next() * SafeGuardMinimum;
        double tmp_weight1 = 0;
        double tmp_weight2 = 0;
        for (int i = 0; i <= half_window_length * 2; ++i) {
            tmp_weight1 += waveform[i];
            tmp_weight2 += window[i];
        }
        double weighting_coefficient = tmp_weight1 / tmp_weight2;
        for (int i = 0; i <= half_window_length * 2; ++i)
            waveform[i] -= window[i] * weighting_coefficient;

        delete[] base_index;
        delete[] window;
    }

    //-----------------------------------------------------------------------------
    // CheapTrickGeneralBody() calculates a spectral envelope at a temporal
    // position. This function is only used in CheapTrick().
    // Caution:
    //   forward_fft is allocated in advance to speed up the processing.
    //-----------------------------------------------------------------------------
    void CheapTrickGeneralBody(const double *x, int fs,
                               double current_f0, int fft_size, double q1,
                               const ForwardRealFFT *forward_real_fft,
                               const InverseRealFFT *inverse_real_fft, double *spectral_envelope,
                               GaussianNoiseGenerator *randn) {
        // F0-adaptive windowing
        GetWindowedWaveform(x, fs, current_f0, forward_real_fft, randn);

        // Calculate power spectrum with DC correction
        // Note: The calculated power spectrum is stored in an array for waveform.
        // In this imprementation, power spectrum is transformed by FFT (NOT IFFT).
        // However, the same result is obtained.
        // This is tricky but important for simple implementation.
        GetPowerSpectrum(fs, current_f0, fft_size, forward_real_fft);

        // Smoothing of the power (linear axis)
        // forward_real_fft.waveform is the power spectrum.
        LinearSmoothing(forward_real_fft->waveform, current_f0 * 2.0 / 3.0, fs, fft_size, forward_real_fft->waveform);

        // Add infinitesimal noise
        // This is a safeguard to avoid including zero in the spectrum.
        AddInfinitesimalNoise(forward_real_fft->waveform, fft_size, forward_real_fft->waveform, randn);

        // Smoothing (log axis) and spectral recovery on the cepstrum domain.
        SmoothingWithRecovery(current_f0, fs, fft_size, q1, forward_real_fft, inverse_real_fft, spectral_envelope);
    }
    double GetF0FloorForCheapTrick(int fs, int fft_size) {
        return 3.0 * fs / (fft_size - 3.0);
    }
}  // namespace

AnalyzePeriodicityWithCheapTrick::AnalyzePeriodicityWithCheapTrick(unsigned int samplingFrequency)
        : fftSize(uzume::dsp::world::fftSize(samplingFrequency)), samplingFrequency(samplingFrequency),
          forwardRealFft({0}), inverseRealFft({0}),
          randn() {
    InitializeForwardRealFFT(fftSize, &forwardRealFft);
    InitializeInverseRealFFT(fftSize, &inverseRealFft);

}

AnalyzePeriodicityWithCheapTrick::~AnalyzePeriodicityWithCheapTrick() noexcept {
    DestroyForwardRealFFT(&forwardRealFft);
}

bool AnalyzePeriodicityWithCheapTrick::operator()(Spectrum *output, const InstantWaveform *input) {
    double f0 = input->f0 <= GetF0FloorForCheapTrick(input->samplingFrequency, output->fftSize) ? DefaultF0 : input->f0;
    CheapTrickGeneralBody(input->data, input->samplingFrequency, f0, output->fftSize,
                          DefaultQ1, &forwardRealFft, &inverseRealFft, output->periodicSpectrum, &randn);
    return true;
}
