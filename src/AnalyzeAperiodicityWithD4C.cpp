// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include "AnalyzeAperiodicityWithD4C.hpp"
#include <algorithm>
#include <cmath>
#include "constant.hpp"
#include "util.hpp"

using namespace uzume::dsp;

namespace {
    //-----------------------------------------------------------------------------
    // SetParametersForGetWindowedWaveform()
    //-----------------------------------------------------------------------------
    void SetParametersForGetWindowedWaveform(int half_window_length,
                                             int x_length, double current_position, int fs, double current_f0,
                                             WindowFunctionType wt, double window_length_ratio, int *base_index,
                                             int *safe_index, double *window) {
        for (int i = -half_window_length; i <= half_window_length; ++i)
            base_index[i + half_window_length] = i;
        int origin = matlab_round(current_position * fs + 0.001);
        for (int i = 0; i <= half_window_length * 2; ++i)
            safe_index[i] = myMin<int>(x_length - 1, myMax<int>(0, origin + base_index[i]));

        // Designing of the window function
        double position;
        switch (wt) {
            case uzume::dsp::Hanning:
                for (int i = 0; i <= half_window_length * 2; ++i) {
                    position = (2.0 * base_index[i] / window_length_ratio) / fs;
                    window[i] = 0.5 * cos(Pi * position * current_f0) + 0.5;
                }
                break;
            case uzume::dsp::Blackman:
            default:
                for (int i = 0; i <= half_window_length * 2; ++i) {
                    position = (2.0 * base_index[i] / window_length_ratio) / fs;
                    window[i] = 0.42 + 0.5 * cos(Pi * position * current_f0) + 0.08 * cos(Pi * position * current_f0 * 2);
                }
        }
    }

    //-----------------------------------------------------------------------------
    // GetWindowedWaveform() windows the waveform by F0-adaptive window
    // In the variable window_type, 1: hanning, 2: blackman
    //-----------------------------------------------------------------------------
    void GetWindowedWaveform(const double *x, int x_length, int fs,
                             double current_f0, double current_position, WindowFunctionType wt,
                             double window_length_ratio, double *waveform, GaussianNoiseGenerator *randn) {
        int half_window_length =
                matlab_round(window_length_ratio * fs / current_f0 / 2.0);

        int *base_index = new int[half_window_length * 2 + 1];
        int *safe_index = new int[half_window_length * 2 + 1];
        double *window  = new double[half_window_length * 2 + 1];

        SetParametersForGetWindowedWaveform(half_window_length, x_length,
                                            current_position, fs, current_f0, wt, window_length_ratio,
                                            base_index, safe_index, window);

        // F0-adaptive windowing
        for (int i = 0; i <= half_window_length * 2; ++i)
            waveform[i] =
                    x[safe_index[i]] * window[i] + randn->next() * SafeGuardMinimum;

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
        delete[] safe_index;
        delete[] window;
    }

    //-----------------------------------------------------------------------------
    // GetCentroid() calculates the energy centroid (see the book, time-frequency
    // analysis written by L. Cohen).
    //-----------------------------------------------------------------------------
    void GetCentroid(const double *x, int x_length, int fs,
                     double current_f0, int fft_size, double current_position,
                     const ForwardRealFFT *forward_real_fft, double *centroid, GaussianNoiseGenerator *randn) {
        for (int i = 0; i < fft_size; ++i) forward_real_fft->waveform[i] = 0.0;
        GetWindowedWaveform(x, x_length, fs, current_f0, current_position, Blackman, 4.0, forward_real_fft->waveform, randn);
        double power = 0.0;
        for (int i = 0; i <= matlab_round(2.0 * fs / current_f0) * 2; ++i)
            power += forward_real_fft->waveform[i] * forward_real_fft->waveform[i];
        for (int i = 0; i <= matlab_round(2.0 * fs / current_f0) * 2; ++i)
            forward_real_fft->waveform[i] /= sqrt(power);

        fft_execute(forward_real_fft->forward_fft);
        double *tmp_real = new double[fft_size / 2 + 1];
        double *tmp_imag = new double[fft_size / 2 + 1];
        for (int i = 0; i <= fft_size / 2; ++i) {
            tmp_real[i] = forward_real_fft->spectrum[i][0];
            tmp_imag[i] = forward_real_fft->spectrum[i][1];
        }

        for (int i = 0; i < fft_size; ++i)
            forward_real_fft->waveform[i] *= i + 1.0;
        fft_execute(forward_real_fft->forward_fft);
        for (int i = 0; i <= fft_size / 2; ++i)
            centroid[i] = forward_real_fft->spectrum[i][0] * tmp_real[i] +
                          tmp_imag[i] * forward_real_fft->spectrum[i][1];

        delete[] tmp_real;
        delete[] tmp_imag;
    }

    //-----------------------------------------------------------------------------
    // GetStaticCentroid() calculates the temporally static energy centroid.
    // Basic idea was proposed by H. Kawahara.
    //-----------------------------------------------------------------------------
    void GetStaticCentroid(const double *x, int x_length, int fs,
                           double current_f0, int fft_size, double current_position,
                           const ForwardRealFFT *forward_real_fft, double *static_centroid, GaussianNoiseGenerator *randn) {
        double *centroid1 = new double[fft_size / 2 + 1];
        double *centroid2 = new double[fft_size / 2 + 1];

        GetCentroid(x, x_length, fs, current_f0, fft_size,
                    current_position - 0.25 / current_f0, forward_real_fft, centroid1, randn);
        GetCentroid(x, x_length, fs, current_f0, fft_size,
                    current_position + 0.25 / current_f0, forward_real_fft, centroid2, randn);

        for (int i = 0; i <= fft_size / 2; ++i)
            static_centroid[i] = centroid1[i] + centroid2[i];

        DCCorrection(static_centroid, current_f0, fs, fft_size, static_centroid);
        delete[] centroid1;
        delete[] centroid2;
    }

//-----------------------------------------------------------------------------
// GetSmoothedPowerSpectrum() calculates the smoothed power spectrum.
// The parameters used for smoothing are optimized in davance.
//-----------------------------------------------------------------------------
    void GetSmoothedPowerSpectrum(const double *x, int x_length, int fs,
                                  double current_f0, int fft_size, double current_position,
                                  const ForwardRealFFT *forward_real_fft, double *smoothed_power_spectrum,
                                  GaussianNoiseGenerator *randn) {
        for (int i = 0; i < fft_size; ++i) forward_real_fft->waveform[i] = 0.0;
        GetWindowedWaveform(x, x_length, fs, current_f0, current_position, Hanning, 4.0, forward_real_fft->waveform, randn);

        fft_execute(forward_real_fft->forward_fft);
        for (int i = 0; i <= fft_size / 2; ++i)
            smoothed_power_spectrum[i] =
                    forward_real_fft->spectrum[i][0] * forward_real_fft->spectrum[i][0] +
                    forward_real_fft->spectrum[i][1] * forward_real_fft->spectrum[i][1];
        DCCorrection(smoothed_power_spectrum, current_f0, fs, fft_size,
                     smoothed_power_spectrum);
        LinearSmoothing(smoothed_power_spectrum, current_f0, fs, fft_size,
                        smoothed_power_spectrum);
    }

//-----------------------------------------------------------------------------
// GetStaticGroupDelay() calculates the temporally static group delay.
// This is the fundamental parameter in D4C.
//-----------------------------------------------------------------------------
    static void GetStaticGroupDelay(const double *static_centroid,
                                    const double *smoothed_power_spectrum, int fs, double f0,
                                    int fft_size, double *static_group_delay) {
        for (int i = 0; i <= fft_size / 2; ++i)
            static_group_delay[i] = static_centroid[i] / smoothed_power_spectrum[i];
        LinearSmoothing(static_group_delay, f0 / 2.0, fs, fft_size, static_group_delay);

        double *smoothed_group_delay = new double[fft_size / 2 + 1];
        LinearSmoothing(static_group_delay, f0, fs, fft_size, smoothed_group_delay);

        for (int i = 0; i <= fft_size / 2; ++i)
            static_group_delay[i] -= smoothed_group_delay[i];

        delete[] smoothed_group_delay;
    }

//-----------------------------------------------------------------------------
// GetCoarseAperiodicity() calculates the aperiodicity in multiples of 3 kHz.
// The upper limit is given based on the sampling frequency.
//-----------------------------------------------------------------------------
    static void GetCoarseAperiodicity(const double *static_group_delay, int fs,
                                      int fft_size, int number_of_aperiodicities, const double *window,
                                      int window_length, const ForwardRealFFT *forward_real_fft,
                                      double *coarse_aperiodicity) {
        int boundary =
                matlab_round(fft_size * 8.0 / window_length);
        int half_window_length = window_length / 2;

        for (int i = 0; i < fft_size; ++i) forward_real_fft->waveform[i] = 0.0;

        double *power_spectrum = new double[fft_size / 2 + 1];
        int center;
        for (int i = 0; i < number_of_aperiodicities; ++i) {
            center = static_cast<int>(FrequencyInterval * (i + 1) * fft_size / fs);
            for (int j = 0; j <= half_window_length * 2; ++j)
                forward_real_fft->waveform[j] = static_group_delay[center - half_window_length + j] * window[j];
            fft_execute(forward_real_fft->forward_fft);
            for (int j = 0 ; j <= fft_size / 2; ++j)
                power_spectrum[j] =
                        forward_real_fft->spectrum[j][0] * forward_real_fft->spectrum[j][0] +
                        forward_real_fft->spectrum[j][1] * forward_real_fft->spectrum[j][1];
            std::sort(power_spectrum, power_spectrum + fft_size / 2 + 1);
            for (int j = 1 ; j <= fft_size / 2; ++j)
                power_spectrum[j] += power_spectrum[j - 1];
            coarse_aperiodicity[i] =
                    10 * log10(power_spectrum[fft_size / 2 - boundary - 1] / power_spectrum[fft_size / 2]);
        }
        delete[] power_spectrum;
    }
    //-----------------------------------------------------------------------------
    // D4CGeneralBody() calculates a spectral envelope at a temporal
    // position. This function is only used in D4C().
    // Caution:
    //   forward_fft is allocated in advance to speed up the processing.
    //-----------------------------------------------------------------------------
    void D4CGeneralBody(const double *x, int x_length, int fs,
                        double current_f0, int fft_size, double current_position,
                        int number_of_aperiodicities, const double *window, int window_length,
                        const ForwardRealFFT *forward_real_fft, double *coarse_aperiodicity,
                        GaussianNoiseGenerator *randn) {
        double *static_centroid = new double[fft_size / 2 + 1];
        double *smoothed_power_spectrum = new double[fft_size / 2 + 1];
        double *static_group_delay = new double[fft_size / 2 + 1];
        GetStaticCentroid(x, x_length, fs, current_f0, fft_size, current_position, forward_real_fft, static_centroid, randn);
        GetSmoothedPowerSpectrum(x, x_length, fs, current_f0, fft_size,
                                 current_position, forward_real_fft, smoothed_power_spectrum, randn);
        GetStaticGroupDelay(static_centroid, smoothed_power_spectrum,
                            fs, current_f0, fft_size, static_group_delay);

        GetCoarseAperiodicity(static_group_delay, fs, fft_size,
                              number_of_aperiodicities, window, window_length, forward_real_fft,
                              coarse_aperiodicity);

        // Revision of the result based on the F0
        for (int i = 0; i < number_of_aperiodicities; ++i)
            coarse_aperiodicity[i] = myMin<double>(0.0,coarse_aperiodicity[i] + (current_f0 - 100) / 50.0);

        delete[] static_centroid;
        delete[] smoothed_power_spectrum;
        delete[] static_group_delay;
    }
    void GetAperiodicity(const double *coarse_frequency_axis,
                         const double *coarse_aperiodicity, int number_of_aperiodicities,
                         const double *frequency_axis, int fft_size, double *aperiodicity) {
        interp1(coarse_frequency_axis, coarse_aperiodicity,
                number_of_aperiodicities + 2, frequency_axis, fft_size / 2 + 1,
                aperiodicity);
        for (int i = 0; i <= fft_size / 2; ++i)
            aperiodicity[i] = pow(10.0, aperiodicity[i] / 20.0);
    }

    void InitializeAperiodicSpectrum(double *x, int x_length) {
        for (int i = 0; i < x_length; i++) {
            x[i] = 1.0 - SafeGuardMinimum;
        }
    }

    const double Boundary0F0 = 100.0;
    const double Boundary1F0 = 4000.0;
    const double Boundary2F0 = 7900.0;

    double D4CLoveTrainSub(const double *x, int fs, int x_length,
                           double current_f0, double current_position, int fft_size,
                           int boundary0, int boundary1, int boundary2,
                           ForwardRealFFT *forward_real_fft, GaussianNoiseGenerator *randn) {
        if(current_f0 == 0.0) {
            return 0.0;
        }
        double *power_spectrum = new double[fft_size];

        int window_length = matlab_round(1.5 * fs / current_f0) * 2 + 1;
        GetWindowedWaveform(x, x_length, fs, current_f0, current_position, Blackman, 3.0, forward_real_fft->waveform, randn);

        for (int i = window_length; i < fft_size; ++i)
            forward_real_fft->waveform[i] = 0.0;
        fft_execute(forward_real_fft->forward_fft);

        for (int i = 0; i <= boundary0; ++i) power_spectrum[i] = 0.0;
        for (int i = boundary0 + 1; i < fft_size / 2 + 1; ++i)
            power_spectrum[i] =
                    forward_real_fft->spectrum[i][0] * forward_real_fft->spectrum[i][0] +
                    forward_real_fft->spectrum[i][1] * forward_real_fft->spectrum[i][1];
        for (int i = boundary0; i <= boundary2; ++i)
            power_spectrum[i] += +power_spectrum[i - 1];

        double aperiodicity0 = power_spectrum[boundary1] / power_spectrum[boundary2];
        delete[] power_spectrum;
        return aperiodicity0;
    }
}  // namespace

AnalyzeAperiodicityWithD4C::AnalyzeAperiodicityWithD4C(unsigned int fftSize, unsigned int samplingFrequency)
        : fftSize(fftSize), samplingFrequency(samplingFrequency),
          forwardRealFft({0}), inverseRealFft({0}),
          coarseAperiodicity(nullptr), coarseFrequencyAxis(nullptr), frequencyAxis(nullptr), nuttallWindow(nullptr),
          forwardRealFFtForD4CLoveTrain({0}),
          boundary0((int)ceil(Boundary0F0 * fftSize / samplingFrequency)),
          boundary1((int)ceil(Boundary1F0 * fftSize / samplingFrequency)),
          boundary2((int)ceil(Boundary2F0 * fftSize / samplingFrequency)),
          randn() {
    int fs4D4C = fftSizeForD4C();
    InitializeForwardRealFFT(fs4D4C, &forwardRealFft);
    InitializeInverseRealFFT(fs4D4C, &inverseRealFft);


    InitializeForwardRealFFT(fftSizeForD4CLoveTrain(), &forwardRealFFtForD4CLoveTrain);
    int numAp = numberOfAperiodicities();

    coarseAperiodicity = new double[numAp + 2];
    coarseAperiodicity[0] = -60.0;
    coarseAperiodicity[numAp + 1] = -SafeGuardMinimum;

    coarseFrequencyAxis = new double[numAp + 2];
    for(int i = 0; i <= numAp; i++) {
        coarseFrequencyAxis[i] = i * FrequencyInterval;
    }
    coarseFrequencyAxis[numAp + 1] = samplingFrequency / 2.0;

    int windowSize = (int)(FrequencyInterval * fs4D4C / samplingFrequency) * 2 + 1;
    nuttallWindow = new double[windowSize];
    NuttallWindow(windowSize, nuttallWindow);

    frequencyAxis = new double[fftSize / 2 + 1];
    for(unsigned int i = 0; i <= fftSize / 2; i++) {
        frequencyAxis[i] = (double)i * samplingFrequency / fftSize;
    }
}

AnalyzeAperiodicityWithD4C::~AnalyzeAperiodicityWithD4C() noexcept {
    DestroyForwardRealFFT(&forwardRealFft);
    DestroyInverseRealFFT(&inverseRealFft);
    DestroyForwardRealFFT(&forwardRealFFtForD4CLoveTrain);
    delete[] coarseAperiodicity;
    delete[] coarseFrequencyAxis;
    delete[] nuttallWindow;
    delete[] frequencyAxis;
}

int AnalyzeAperiodicityWithD4C::numberOfAperiodicities() const {
    return (int)(myMin<double>(UpperLimit, samplingFrequency / 2.0 - FrequencyInterval) / FrequencyInterval);;
}

int AnalyzeAperiodicityWithD4C::fftSizeForD4C() const {
    return (int)(pow(2.0, 1.0 + (int)(log(4.0 * samplingFrequency / FloorF0D4C + 1) / Log2)));
}

int AnalyzeAperiodicityWithD4C::fftSizeForD4CLoveTrain() const {
    return (int)(pow(2.0, 1.0 + (int)(log(3.0 * samplingFrequency / LowestF0D4CLoveTrain + 1) / Log2)));
}

int AnalyzeAperiodicityWithD4C::nuttallWindowSize() const {
    return (int)(FrequencyInterval * fftSizeForD4C() / samplingFrequency) * 2 + 1;
}

bool AnalyzeAperiodicityWithD4C::operator()(Spectrum *output, const InstantWaveform *input) {
    InitializeAperiodicSpectrum(output->aperiodicSpectrum, output->fftSize);
    int fftSizeForAperiodicity0 = fftSizeForD4CLoveTrain();
    double aperiodicity0 = D4CLoveTrainSub(input->data, samplingFrequency, input->length, input->f0, 0.0, fftSizeForAperiodicity0,
                                           Boundary0F0 * fftSizeForAperiodicity0 / samplingFrequency,
                                           Boundary1F0 * fftSizeForAperiodicity0 / samplingFrequency,
                                           Boundary2F0 * fftSizeForAperiodicity0 / samplingFrequency,
                                           &forwardRealFFtForD4CLoveTrain, &randn);
    if (input->f0 == 0 || aperiodicity0 <= Threshold) {
        return true;
    }
    D4CGeneralBody(input->data, input->length, samplingFrequency, myMax<double>(FloorF0D4C, input->f0),
                   fftSizeForD4C(), 0.0, numberOfAperiodicities(), nuttallWindow,
                   nuttallWindowSize(), &forwardRealFft, &coarseAperiodicity[1], &randn);

    // Linear interpolation to convert the coarse aperiodicity into its
    // spectral representation.
    GetAperiodicity(coarseFrequencyAxis, coarseAperiodicity, numberOfAperiodicities(), frequencyAxis, fftSize, output->aperiodicSpectrum);
    return true;
}
