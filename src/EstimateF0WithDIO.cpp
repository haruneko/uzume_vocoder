// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include <cmath>
#include "constant.hpp"
#include "EstimateF0WithDIO.hpp"
#include "fft.hpp"
#include "util.hpp"

using namespace uzume::dsp;

//-----------------------------------------------------------------------------
// struct for GetFourZeroCrossingIntervals()
// "negative" means "zero-crossing point going from positive to negative"
// "positive" means "zero-crossing point going from negative to positive"
//-----------------------------------------------------------------------------
typedef struct {
    double *negative_interval_locations;
    double *negative_intervals;
    int number_of_negatives;
    double *positive_interval_locations;
    double *positive_intervals;
    int number_of_positives;
    double *peak_interval_locations;
    double *peak_intervals;
    int number_of_peaks;
    double *dip_interval_locations;
    double *dip_intervals;
    int number_of_dips;
} ZeroCrossings;

// From DIO
namespace {

    int GetSuitableFFTSize(int sample) {
        return (int) (pow(2.0, (int) (log(static_cast<double>(sample)) / Log2) + 1.0));
    }
    int GetSamplesForDIO(int fs, int x_length, double frame_period) {
        return static_cast<int>(1000.0 * x_length / fs / frame_period) + 1;
    }
//-----------------------------------------------------------------------------
// FilterForDecimate() calculates the coefficients of low-pass filter and
// carries out the filtering. This function is only used for decimate().
//-----------------------------------------------------------------------------
    static void FilterForDecimate(const double *x, int x_length, int r, double *y) {
        double a[3], b[2];  // filter Coefficients
        switch (r) {
            case 11:  // fs : 44100 (default)
                a[0] = 2.450743295230728;
                a[1] = -2.06794904601978;
                a[2] = 0.59574774438332101;
                b[0] = 0.0026822508007163792;
                b[1] = 0.0080467524021491377;
                break;
            case 12:  // fs : 48000
                a[0] = 2.4981398605924205;
                a[1] = -2.1368928194784025;
                a[2] = 0.62187513816221485;
                b[0] = 0.0021097275904709001;
                b[1] = 0.0063291827714127002;
                break;
            case 10:
                a[0] = 2.3936475118069387;
                a[1] = -1.9873904075111861;
                a[2] = 0.5658879979027055;
                b[0] = 0.0034818622251927556;
                b[1] = 0.010445586675578267;
                break;
            case 9:
                a[0] = 2.3236003491759578;
                a[1] = -1.8921545617463598;
                a[2] = 0.53148928133729068;
                b[0] = 0.0046331164041389372;
                b[1] = 0.013899349212416812;
                break;
            case 8:  // fs : 32000
                a[0] = 2.2357462340187593;
                a[1] = -1.7780899984041358;
                a[2] = 0.49152555365968692;
                b[0] = 0.0063522763407111993;
                b[1] = 0.019056829022133598;
                break;
            case 7:
                a[0] = 2.1225239019534703;
                a[1] = -1.6395144861046302;
                a[2] = 0.44469707800587366;
                b[0] = 0.0090366882681608418;
                b[1] = 0.027110064804482525;
                break;
            case 6:  // fs : 24000 and 22050
                a[0] = 1.9715352749512141;
                a[1] = -1.4686795689225347;
                a[2] = 0.3893908434965701;
                b[0] = 0.013469181309343825;
                b[1] = 0.040407543928031475;
                break;
            case 5:
                a[0] = 1.7610939654280557;
                a[1] = -1.2554914843859768;
                a[2] = 0.3237186507788215;
                b[0] = 0.021334858522387423;
                b[1] = 0.06400457556716227;
                break;
            case 4:  // fs : 16000
                a[0] = 1.4499664446880227;
                a[1] = -0.98943497080950582;
                a[2] = 0.24578252340690215;
                b[0] = 0.036710750339322612;
                b[1] = 0.11013225101796784;
                break;
            case 3:
                a[0] = 0.95039378983237421;
                a[1] = -0.67429146741526791;
                a[2] = 0.15412211621346475;
                b[0] = 0.071221945171178636;
                b[1] = 0.21366583551353591;
                break;
            case 2:  // fs : 8000
                a[0] = 0.041156734567757189;
                a[1] = -0.42599112459189636;
                a[2] = 0.041037215479961225;
                b[0] = 0.16797464681802227;
                b[1] = 0.50392394045406674;
                break;
            default:
                a[0] = 0.0;
                a[1] = 0.0;
                a[2] = 0.0;
                b[0] = 0.0;
                b[1] = 0.0;
        }

        // Filtering on time domain.
        double w[3] = {0.0, 0.0, 0.0};
        double wt;
        for (int i = 0; i < x_length; ++i) {
            wt = x[i] + a[0] * w[0] + a[1] * w[1] + a[2] * w[2];
            y[i] = b[0] * wt + b[1] * w[0] + b[1] * w[1] + b[0] * w[2];
            w[2] = w[1];
            w[1] = w[0];
            w[0] = wt;
        }
    }

    void decimate(const double *x, int x_length, int r, double *y) {
        const int kNFact = 9;
        double *tmp1 = new double[x_length + kNFact * 2];
        double *tmp2 = new double[x_length + kNFact * 2];

        for (int i = 0; i < kNFact; ++i) tmp1[i] = 2 * x[0] - x[kNFact - i];
        for (int i = kNFact; i < kNFact + x_length; ++i) tmp1[i] = x[i - kNFact];
        for (int i = kNFact + x_length; i < 2 * kNFact + x_length; ++i)
            tmp1[i] = 2 * x[x_length - 1] - x[x_length - 2 - (i - (kNFact + x_length))];

        FilterForDecimate(tmp1, 2 * kNFact + x_length, r, tmp2);
        for (int i = 0; i < 2 * kNFact + x_length; ++i)
            tmp1[i] = tmp2[2 * kNFact + x_length - i - 1];
        FilterForDecimate(tmp1, 2 * kNFact + x_length, r, tmp2);
        for (int i = 0; i < 2 * kNFact + x_length; ++i)
            tmp1[i] = tmp2[2 * kNFact + x_length - i - 1];

        int nout = (x_length - 1) / r + 1;
        int nbeg = r - r * nout + x_length;

        int count = 0;
        for (int i = nbeg; i < x_length + kNFact; i += r)
            y[count++] = tmp1[i + kNFact - 1];

        delete[] tmp1;
        delete[] tmp2;
    }

//-----------------------------------------------------------------------------
// DesignLowCutFilter() calculates the coefficients the filter.
//-----------------------------------------------------------------------------
    static void DesignLowCutFilter(int N, int fft_size, double *low_cut_filter) {
        for (int i = 1; i <= N; ++i)
            low_cut_filter[i - 1] = 0.5 - 0.5 * cos(i * 2.0 * Pi / (N + 1));
        for (int i = N; i < fft_size; ++i) low_cut_filter[i] = 0.0;
        double sum_of_amplitude = 0.0;
        for (int i = 0; i < N; ++i) sum_of_amplitude += low_cut_filter[i];
        for (int i = 0; i < N; ++i)
            low_cut_filter[i] = -low_cut_filter[i] / sum_of_amplitude;
        for (int i = 0; i < (N - 1) / 2; ++i)
            low_cut_filter[fft_size - (N - 1) / 2 + i] = low_cut_filter[i];
        for (int i = 0; i < N; ++i)
            low_cut_filter[i] = low_cut_filter[i + (N - 1) / 2];
        low_cut_filter[0] += 1.0;
    }

//-----------------------------------------------------------------------------
// GetSpectrumForEstimation() calculates the spectrum for estimation.
// This function carries out downsampling to speed up the estimation process
// and calculates the spectrum of the downsampled signal.
//-----------------------------------------------------------------------------
    static void GetSpectrumForEstimation(const double *x, int x_length,
                                         int y_length, double actual_fs, int fft_size, int decimation_ratio,
                                         fft_complex *y_spectrum) {
        double *y = new double[fft_size];

        // Initialization
        for (int i = 0; i < fft_size; ++i) y[i] = 0.0;

        // Downsampling
        if (decimation_ratio != 1)
            decimate(x, x_length, decimation_ratio, y);
        else
            for (int i = 0; i < x_length; ++i) y[i] = x[i];

        // Removal of the DC component (y = y - mean value of y)
        double mean_y = 0.0;
        for (int i = 0; i < y_length; ++i) mean_y += y[i];
        mean_y /= y_length;
        for (int i = 0; i < y_length; ++i) y[i] -= mean_y;
        for (int i = y_length; i < fft_size; ++i) y[i] = 0.0;

        fft_plan forwardFFT = fft_plan_dft_r2c_1d(fft_size, y, y_spectrum, FFT_ESTIMATE);
        fft_execute(forwardFFT);

        // Low cut filtering (from 0.1.4). Cut off frequency is 50.0 Hz.
        int cutoff_in_sample = matlab_round(actual_fs / CutOff);
        DesignLowCutFilter(cutoff_in_sample * 2 + 1, fft_size, y);

        fft_complex *filter_spectrum = new fft_complex[fft_size];
        forwardFFT.c_out = filter_spectrum;
        fft_execute(forwardFFT);

        double tmp = 0;
        for (int i = 0; i <= fft_size / 2; ++i) {
            // Complex number multiplications.
            tmp = y_spectrum[i][0] * filter_spectrum[i][0] - y_spectrum[i][1] * filter_spectrum[i][1];
            y_spectrum[i][1] = y_spectrum[i][0] * filter_spectrum[i][1] + y_spectrum[i][1] * filter_spectrum[i][0];
            y_spectrum[i][0] = tmp;
        }

        fft_destroy_plan(forwardFFT);
        delete[] y;
        delete[] filter_spectrum;
    }

//-----------------------------------------------------------------------------
// GetBestF0Contour() calculates the best data contour based on scores of
// all candidates. The F0 with highest score is selected.
//-----------------------------------------------------------------------------
    static void GetBestF0Contour(int f0_length,
                                 const double *const *f0_candidates, const double *const *f0_scores,
                                 int number_of_bands, double *best_f0_contour) {
        double tmp;
        for (int i = 0; i < f0_length; ++i) {
            tmp = f0_scores[0][i];
            best_f0_contour[i] = f0_candidates[0][i];
            for (int j = 1; j < number_of_bands; ++j) {
                if (tmp > f0_scores[j][i]) {
                    tmp = f0_scores[j][i];
                    best_f0_contour[i] = f0_candidates[j][i];
                }
            }
        }
    }

//-----------------------------------------------------------------------------
// FixStep1() is the 1st step of the postprocessing.
// This function eliminates the unnatural change of data based on allowed_range.
//-----------------------------------------------------------------------------
    static void FixStep1(const double *best_f0_contour, int f0_length,
                         int voice_range_minimum, double allowed_range, double *f0_step1) {
        double *f0_base = new double[f0_length];
        // Initialization
        for (int i = 0; i < voice_range_minimum; ++i) f0_base[i] = 0.0;
        for (int i = voice_range_minimum; i < f0_length - voice_range_minimum; ++i)
            f0_base[i] = best_f0_contour[i];
        for (int i = f0_length - voice_range_minimum; i < f0_length; ++i)
            f0_base[i] = 0.0;

        // Processing to prevent the jumping of data
        for (int i = 0; i < voice_range_minimum; ++i) f0_step1[i] = 0.0;
        for (int i = voice_range_minimum; i < f0_length; ++i)
            f0_step1[i] =
                    fabs((f0_base[i] - f0_base[i - 1]) / (SafeGuardMinimum + f0_base[i])) < allowed_range ? f0_base[i] : 0.0;

        delete[] f0_base;
    }

//-----------------------------------------------------------------------------
// FixStep2() is the 2nd step of the postprocessing.
// This function eliminates the suspected data in the anlaut and auslaut.
//-----------------------------------------------------------------------------
    static void FixStep2(const double *f0_step1, int f0_length, int voice_range_minimum, double *f0_step2) {
        for (int i = 0; i < f0_length; ++i) f0_step2[i] = f0_step1[i];

        int center = (voice_range_minimum - 1) / 2;
        for (int i = center; i < f0_length - center; ++i) {
            for (int j = -center; j <= center; ++j) {
                if (f0_step1[i + j] == 0) {
                    f0_step2[i] = 0.0;
                    break;
                }
            }
        }
    }

//-----------------------------------------------------------------------------
// GetNumberOfVoicedSections() counts the number of voiced sections.
//-----------------------------------------------------------------------------
    static void GetNumberOfVoicedSections(const double *f0, int f0_length,
                                          int *positive_index, int *negative_index, int *positive_count,
                                          int *negative_count) {
        *positive_count = *negative_count = 0;
        for (int i = 1; i < f0_length; ++i)
            if (f0[i] == 0 && f0[i - 1] != 0)
                negative_index[(*negative_count)++] = i - 1;
            else if (f0[i - 1] == 0 && f0[i] != 0)
                positive_index[(*positive_count)++] = i;
    }

//-----------------------------------------------------------------------------
// SelectOneF0() corrects the data[current_index] based on
// data[current_index + sign].
//-----------------------------------------------------------------------------
    static double SelectBestF0(double current_f0, double past_f0,
                               const double *const *f0_candidates, int number_of_candidates,
                               int target_index, double allowed_range) {
        double reference_f0 = (current_f0 * 3.0 - past_f0) / 2.0;

        double minimum_error = fabs(reference_f0 - f0_candidates[0][target_index]);
        double best_f0 = f0_candidates[0][target_index];

        double current_error;
        for (int i = 1; i < number_of_candidates; ++i) {
            current_error = fabs(reference_f0 - f0_candidates[i][target_index]);
            if (current_error < minimum_error) {
                minimum_error = current_error;
                best_f0 = f0_candidates[i][target_index];
            }
        }
        if (fabs(1.0 - best_f0 / reference_f0) > allowed_range)
            return 0.0;
        return best_f0;
    }

//-----------------------------------------------------------------------------
// FixStep3() is the 3rd step of the postprocessing.
// This function corrects the data candidates from backward to forward.
//-----------------------------------------------------------------------------
    static void FixStep3(const double *f0_step2, int f0_length,
                         const double *const *f0_candidates, int number_of_candidates,
                         double allowed_range, const int *negative_index, int negative_count,
                         double *f0_step3) {
        for (int i = 0; i < f0_length; i++) f0_step3[i] = f0_step2[i];

        int limit;
        for (int i = 0; i < negative_count; ++i) {
            limit = i == negative_count - 1 ? f0_length - 1 : negative_index[i + 1];
            for (int j = negative_index[i]; j < limit; ++j) {
                f0_step3[j + 1] =
                        SelectBestF0(f0_step3[j], f0_step3[j - 1], f0_candidates,
                                     number_of_candidates, j + 1, allowed_range);
                if (f0_step3[j + 1] == 0) break;
            }
        }
    }

//-----------------------------------------------------------------------------
// FixStep4() is the 4th step of the postprocessing.
// This function corrects the data candidates from forward to backward.
//-----------------------------------------------------------------------------
    static void FixStep4(const double *f0_step3, int f0_length,
                         const double *const *f0_candidates, int number_of_candidates,
                         double allowed_range, const int *positive_index, int positive_count,
                         double *f0_step4) {
        for (int i = 0; i < f0_length; ++i) f0_step4[i] = f0_step3[i];

        int limit;
        for (int i = positive_count - 1; i >= 0; --i) {
            limit = i == 0 ? 1 : positive_index[i - 1];
            for (int j = positive_index[i]; j > limit; --j) {
                f0_step4[j - 1] =
                        SelectBestF0(f0_step4[j], f0_step4[j + 1], f0_candidates,
                                     number_of_candidates, j - 1, allowed_range);
                if (f0_step4[j - 1] == 0) break;
            }
        }
    }

//-----------------------------------------------------------------------------
// FixF0Contour() calculates the definitive data contour based on all data
// candidates. There are four steps.
//-----------------------------------------------------------------------------
    static void FixF0Contour(double frame_period, int number_of_candidates,
                             int fs, const double *const *f0_candidates,
                             const double *best_f0_contour, int f0_length, double f0_floor,
                             double allowed_range, double *fixed_f0_contour) {
        int voice_range_minimum = (int)(0.5 + 1000.0 / frame_period / f0_floor) * 2 + 1;

        if (f0_length <= voice_range_minimum) return;

        double *f0_tmp1 = new double[f0_length];
        double *f0_tmp2 = new double[f0_length];

        FixStep1(best_f0_contour, f0_length, voice_range_minimum, allowed_range, f0_tmp1);
        FixStep2(f0_tmp1, f0_length, voice_range_minimum, f0_tmp2);

        int positive_count, negative_count;
        int *positive_index = new int[f0_length];
        int *negative_index = new int[f0_length];
        GetNumberOfVoicedSections(f0_tmp2, f0_length, positive_index, negative_index, &positive_count, &negative_count);
        FixStep3(f0_tmp2, f0_length, f0_candidates, number_of_candidates,
                 allowed_range, negative_index, negative_count, f0_tmp1);
        FixStep4(f0_tmp1, f0_length, f0_candidates, number_of_candidates,
                 allowed_range, positive_index, positive_count, fixed_f0_contour);

        delete[] f0_tmp1;
        delete[] f0_tmp2;
        delete[] positive_index;
        delete[] negative_index;
    }

//-----------------------------------------------------------------------------
// GetFilteredSignal() calculates the signal that is the convolution of the
// input signal and low-pass filter.
// This function is only used in RawEventByDio()
//-----------------------------------------------------------------------------
    static void GetFilteredSignal(int half_average_length, int fft_size,
                                  const fft_complex *y_spectrum, int y_length, double *filtered_signal) {
        double *low_pass_filter = new double[fft_size];
        // Nuttall window is used as a low-pass filter.
        // Cutoff frequency depends on the window length.
        NuttallWindow(half_average_length * 4, low_pass_filter);
        for (int i = half_average_length * 4; i < fft_size; ++i)
            low_pass_filter[i] = 0.0;

        fft_complex *low_pass_filter_spectrum = new fft_complex[fft_size];
        fft_plan forwardFFT = fft_plan_dft_r2c_1d(fft_size, low_pass_filter, low_pass_filter_spectrum, FFT_ESTIMATE);
        fft_execute(forwardFFT);

        // Convolution
        double tmp = y_spectrum[0][0] * low_pass_filter_spectrum[0][0] - y_spectrum[0][1] * low_pass_filter_spectrum[0][1];
        low_pass_filter_spectrum[0][1] =
                y_spectrum[0][0] * low_pass_filter_spectrum[0][1] + y_spectrum[0][1] * low_pass_filter_spectrum[0][0];
        low_pass_filter_spectrum[0][0] = tmp;
        for (int i = 1; i <= fft_size / 2; ++i) {
            tmp = y_spectrum[i][0] * low_pass_filter_spectrum[i][0] - y_spectrum[i][1] * low_pass_filter_spectrum[i][1];
            low_pass_filter_spectrum[i][1] =
                    y_spectrum[i][0] * low_pass_filter_spectrum[i][1] + y_spectrum[i][1] * low_pass_filter_spectrum[i][0];
            low_pass_filter_spectrum[i][0] = tmp;
            low_pass_filter_spectrum[fft_size - i - 1][0] = low_pass_filter_spectrum[i][0];
            low_pass_filter_spectrum[fft_size - i - 1][1] = low_pass_filter_spectrum[i][1];
        }

        fft_plan inverseFFT = fft_plan_dft_c2r_1d(fft_size, low_pass_filter_spectrum, filtered_signal, FFT_ESTIMATE);
        fft_execute(inverseFFT);

        // Compensation of the delay.
        int index_bias = half_average_length * 2;
        for (int i = 0; i < y_length; ++i)
            filtered_signal[i] = filtered_signal[i + index_bias];

        fft_destroy_plan(inverseFFT);
        fft_destroy_plan(forwardFFT);
        delete[] low_pass_filter_spectrum;
        delete[] low_pass_filter;
    }

//-----------------------------------------------------------------------------
// CheckEvent() returns 1, provided that the input value is over 1.
// This function is for RawEventByDio().
//-----------------------------------------------------------------------------
    static inline int CheckEvent(int x) {
        return x > 0 ? 1 : 0;
    }

//-----------------------------------------------------------------------------
// ZeroCrossingEngine() calculates the zero crossing points from positive to
// negative. Thanks to Custom.Maid http://custom-made.seesaa.net/ (2012/8/19)
//-----------------------------------------------------------------------------
    static int ZeroCrossingEngine(const double *filtered_signal, int y_length,
                                  double fs, double *interval_locations, double *intervals) {
        int *negative_going_points = new int[y_length];

        for (int i = 0; i < y_length - 1; ++i)
            negative_going_points[i] = 0.0 < filtered_signal[i] && filtered_signal[i + 1] <= 0.0 ? i + 1 : 0;
        negative_going_points[y_length - 1] = 0;

        int *edges = new int[y_length];
        int count = 0;
        for (int i = 0; i < y_length; ++i)
            if (negative_going_points[i] > 0)
                edges[count++] = negative_going_points[i];

        if (count < 2) {
            delete[] edges;
            delete[] negative_going_points;
            return 0;
        }

        double *fine_edges = new double[count];
        for (int i = 0; i < count; ++i)
            fine_edges[i] =
                    edges[i] - filtered_signal[edges[i] - 1] /
                               (filtered_signal[edges[i]] - filtered_signal[edges[i] - 1]);

        for (int i = 0; i < count - 1; ++i) {
            intervals[i] = fs / (fine_edges[i + 1] - fine_edges[i]);
            interval_locations[i] = (fine_edges[i] + fine_edges[i + 1]) / 2.0 / fs;
        }

        delete[] fine_edges;
        delete[] edges;
        delete[] negative_going_points;
        return count - 1;
    }

//-----------------------------------------------------------------------------
// GetFourZeroCrossingIntervals() calculates four zero-crossing intervals.
// (1) Zero-crossing going from negative to positive.
// (2) Zero-crossing going from positive to negative.
// (3) Peak, and (4) dip. (3) and (4) are calculated from the zero-crossings of
// the differential of waveform.
//-----------------------------------------------------------------------------
    static void GetFourZeroCrossingIntervals(double *filtered_signal, int y_length, double actual_fs, ZeroCrossings *zero_crossings) {
        // x_length / 4 (old version) is fixed at 2013/07/14
        const int kMaximumNumber = y_length;
        zero_crossings->negative_interval_locations = new double[kMaximumNumber];
        zero_crossings->positive_interval_locations = new double[kMaximumNumber];
        zero_crossings->peak_interval_locations = new double[kMaximumNumber];
        zero_crossings->dip_interval_locations = new double[kMaximumNumber];
        zero_crossings->negative_intervals = new double[kMaximumNumber];
        zero_crossings->positive_intervals = new double[kMaximumNumber];
        zero_crossings->peak_intervals = new double[kMaximumNumber];
        zero_crossings->dip_intervals = new double[kMaximumNumber];

        zero_crossings->number_of_negatives = ZeroCrossingEngine(filtered_signal,
                                                                 y_length, actual_fs,
                                                                 zero_crossings->negative_interval_locations,
                                                                 zero_crossings->negative_intervals);

        for (int i = 0; i < y_length; ++i) filtered_signal[i] = -filtered_signal[i];
        zero_crossings->number_of_positives = ZeroCrossingEngine(filtered_signal,
                                                                 y_length, actual_fs,
                                                                 zero_crossings->positive_interval_locations,
                                                                 zero_crossings->positive_intervals);

        for (int i = 0; i < y_length - 1; ++i)
            filtered_signal[i] =
                    filtered_signal[i] - filtered_signal[i + 1];
        zero_crossings->number_of_peaks = ZeroCrossingEngine(filtered_signal,
                                                             y_length - 1, actual_fs,
                                                             zero_crossings->peak_interval_locations,
                                                             zero_crossings->peak_intervals);

        for (int i = 0; i < y_length - 1; ++i)
            filtered_signal[i] = -filtered_signal[i];
        zero_crossings->number_of_dips = ZeroCrossingEngine(filtered_signal,
                                                            y_length - 1, actual_fs,
                                                            zero_crossings->dip_interval_locations,
                                                            zero_crossings->dip_intervals);
    }

//-----------------------------------------------------------------------------
// GetF0CandidateContourSub() calculates the data candidates and deviations.
// This is the sub-function of GetF0Candidates() and assumes the calculation.
//-----------------------------------------------------------------------------
    static void GetF0CandidateContourSub(
            const double *const *interpolated_f0_set, int f0_length, double f0_floor,
            double f0_ceil, double boundary_f0, double *f0_candidate,
            double *f0_score) {
        for (int i = 0; i < f0_length; ++i) {
            f0_candidate[i] = (interpolated_f0_set[0][i] +
                               interpolated_f0_set[1][i] + interpolated_f0_set[2][i] +
                               interpolated_f0_set[3][i]) / 4.0;

            f0_score[i] = sqrt(((interpolated_f0_set[0][i] - f0_candidate[i]) *
                                (interpolated_f0_set[0][i] - f0_candidate[i]) +
                                (interpolated_f0_set[1][i] - f0_candidate[i]) *
                                (interpolated_f0_set[1][i] - f0_candidate[i]) +
                                (interpolated_f0_set[2][i] - f0_candidate[i]) *
                                (interpolated_f0_set[2][i] - f0_candidate[i]) +
                                (interpolated_f0_set[3][i] - f0_candidate[i]) *
                                (interpolated_f0_set[3][i] - f0_candidate[i])) / 3.0);

            if (f0_candidate[i] > boundary_f0 || f0_candidate[i] < boundary_f0 / 2.0 ||
                f0_candidate[i] > f0_ceil || f0_candidate[i] < f0_floor) {
                f0_candidate[i] = 0.0;
                f0_score[i] = MaximumValue;
            }
        }
    }

//-----------------------------------------------------------------------------
// GetF0CandidateContour() calculates the F0 candidates based on the
// zero-crossings.
//-----------------------------------------------------------------------------
    static void GetF0CandidateContour(const ZeroCrossings *zero_crossings,
                                      double boundary_f0, double f0_floor, double f0_ceil,
                                      const double *temporal_positions, int f0_length,
                                      double *f0_candidate, double *f0_score) {
        if (0 == CheckEvent(zero_crossings->number_of_negatives - 2) *
                 CheckEvent(zero_crossings->number_of_positives - 2) *
                 CheckEvent(zero_crossings->number_of_peaks - 2) *
                 CheckEvent(zero_crossings->number_of_dips - 2)) {
            for (int i = 0; i < f0_length; ++i) {
                f0_score[i] = MaximumValue;
                f0_candidate[i] = 0.0;
            }
            return;
        }

        double *interpolated_f0_set[4];
        for (int i = 0; i < 4; ++i)
            interpolated_f0_set[i] = new double[f0_length];

        interp1(zero_crossings->negative_interval_locations,
                zero_crossings->negative_intervals,
                zero_crossings->number_of_negatives,
                temporal_positions, f0_length, interpolated_f0_set[0]);
        interp1(zero_crossings->positive_interval_locations,
                zero_crossings->positive_intervals,
                zero_crossings->number_of_positives,
                temporal_positions, f0_length, interpolated_f0_set[1]);
        interp1(zero_crossings->peak_interval_locations,
                zero_crossings->peak_intervals, zero_crossings->number_of_peaks,
                temporal_positions, f0_length, interpolated_f0_set[2]);
        interp1(zero_crossings->dip_interval_locations,
                zero_crossings->dip_intervals, zero_crossings->number_of_dips,
                temporal_positions, f0_length, interpolated_f0_set[3]);

        GetF0CandidateContourSub(interpolated_f0_set, f0_length, f0_floor, f0_ceil, boundary_f0, f0_candidate, f0_score);
        for (int i = 0; i < 4; ++i) delete[] interpolated_f0_set[i];
    }

//-----------------------------------------------------------------------------
// DestroyZeroCrossings() frees the memory of array in the struct
//-----------------------------------------------------------------------------
    static void DestroyZeroCrossings(ZeroCrossings *zero_crossings) {
        delete[] zero_crossings->negative_interval_locations;
        delete[] zero_crossings->positive_interval_locations;
        delete[] zero_crossings->peak_interval_locations;
        delete[] zero_crossings->dip_interval_locations;
        delete[] zero_crossings->negative_intervals;
        delete[] zero_crossings->positive_intervals;
        delete[] zero_crossings->peak_intervals;
        delete[] zero_crossings->dip_intervals;
    }

//-----------------------------------------------------------------------------
// GetF0CandidateFromRawEvent() calculates F0 candidate contour in 1-ch signal
//-----------------------------------------------------------------------------
    static void GetF0CandidateFromRawEvent(double boundary_f0, double fs,
                                           const fft_complex *y_spectrum, int y_length, int fft_size, double f0_floor,
                                           double f0_ceil, const double *temporal_positions, int f0_length,
                                           double *f0_score, double *f0_candidate) {
        double *filtered_signal = new double[fft_size];
        GetFilteredSignal(matlab_round(fs / boundary_f0 / 2.0), fft_size, y_spectrum, y_length, filtered_signal);

        ZeroCrossings zero_crossings = {0};
        GetFourZeroCrossingIntervals(filtered_signal, y_length, fs, &zero_crossings);

        GetF0CandidateContour(&zero_crossings, boundary_f0, f0_floor, f0_ceil, temporal_positions, f0_length, f0_candidate, f0_score);

        DestroyZeroCrossings(&zero_crossings);
        delete[] filtered_signal;
    }

//-----------------------------------------------------------------------------
// GetF0CandidatesAndScores() calculates all data candidates and their scores.
//-----------------------------------------------------------------------------
    static void GetF0CandidatesAndScores(const double *boundary_f0_list,
                                         int number_of_bands, double actual_fs, int y_length,
                                         const double *temporal_positions, int f0_length,
                                         const fft_complex *y_spectrum, int fft_size, double f0_floor,
                                         double f0_ceil, double **raw_f0_candidates, double **raw_f0_scores) {
        double *f0_candidate = new double[f0_length];
        double *f0_score = new double[f0_length];

        // Calculation of the acoustics events (zero-crossing)
        for (int i = 0; i < number_of_bands; ++i) {
            GetF0CandidateFromRawEvent(boundary_f0_list[i], actual_fs, y_spectrum,
                                       y_length, fft_size, f0_floor, f0_ceil, temporal_positions, f0_length,
                                       f0_score, f0_candidate);
            for (int j = 0; j < f0_length; ++j) {
                // A way to avoid zero division
                raw_f0_scores[i][j] = f0_score[j] / (f0_candidate[j] + SafeGuardMinimum);
                raw_f0_candidates[i][j] = f0_candidate[j];
            }
        }

        delete[] f0_candidate;
        delete[] f0_score;
    }

//-----------------------------------------------------------------------------
// DioGeneralBody() estimates the F0 based on Distributed Inline-filter
// Operation.
//-----------------------------------------------------------------------------
    static void DioGeneralBody(const double *x, int x_length, int fs,
                               double frame_period, double f0_floor, double f0_ceil,
                               double channels_in_octave, int speed, double allowed_range,
                               double *temporal_positions, double *f0) {
        int number_of_bands = 1 + static_cast<int>(log(f0_ceil / f0_floor) / Log2 * channels_in_octave);
        double *boundary_f0_list = new double[number_of_bands];
        for (int i = 0; i < number_of_bands; ++i)
            boundary_f0_list[i] = f0_floor * pow(2.0, (i + 1) / channels_in_octave);

        // normalization
        int decimation_ratio = myMax<int>(myMin<int>(speed, 12), 1);
        int y_length = (1 + static_cast<int>(x_length / decimation_ratio));
        double actual_fs = static_cast<double>(fs) / decimation_ratio;
        int fft_size = GetSuitableFFTSize(y_length +
                                          matlab_round(actual_fs / CutOff) * 2 + 1 +
                                          (4 * static_cast<int>(1.0 + actual_fs / boundary_f0_list[0] / 2.0)));

        // Calculation of the spectrum used for the data estimation
        fft_complex *y_spectrum = new fft_complex[fft_size];
        GetSpectrumForEstimation(x, x_length, y_length, actual_fs, fft_size, decimation_ratio, y_spectrum);

        double **f0_candidates = new double *[number_of_bands];
        double **f0_scores = new double *[number_of_bands];
        int f0_length = GetSamplesForDIO(fs, x_length, frame_period);
        for (int i = 0; i < number_of_bands; ++i) {
            f0_candidates[i] = new double[f0_length];
            f0_scores[i] = new double[f0_length];
        }

        for (int i = 0; i < f0_length; ++i)
            temporal_positions[i] = i * frame_period / 1000.0;

        GetF0CandidatesAndScores(boundary_f0_list, number_of_bands,
                                 actual_fs, y_length, temporal_positions, f0_length, y_spectrum,
                                 fft_size, f0_floor, f0_ceil, f0_candidates, f0_scores);

        // Selection of the best value based on fundamental-ness.
        // This function is related with SortCandidates() in MATLAB.
        double *best_f0_contour = new double[f0_length];
        GetBestF0Contour(f0_length, f0_candidates, f0_scores, number_of_bands, best_f0_contour);

        // Postprocessing to find the best data-contour.
        FixF0Contour(frame_period, number_of_bands, fs, f0_candidates, best_f0_contour, f0_length, f0_floor, allowed_range, f0);

        delete[] best_f0_contour;
        delete[] y_spectrum;
        for (int i = 0; i < number_of_bands; ++i) {
            delete[] f0_scores[i];
            delete[] f0_candidates[i];
        }
        delete[] f0_scores;
        delete[] f0_candidates;
        delete[] boundary_f0_list;
    }
}

// From StoneMask
namespace {
//-----------------------------------------------------------------------------
// GetBaseIndex() calculates the temporal positions for windowing.
// Since the result includes negative value and the value that exceeds the
// length of the input signal, it must be modified appropriately.
//-----------------------------------------------------------------------------
    static void GetBaseIndex(double current_position, const double *base_time,
                             int base_time_length, int fs, int *index_raw) {
        for (int i = 0; i < base_time_length; ++i)
            index_raw[i] = matlab_round((current_position + base_time[i]) * fs);
    }

//-----------------------------------------------------------------------------
// GetMainWindow() generates the window function.
//-----------------------------------------------------------------------------
    static void GetMainWindow(double current_position, const int *index_raw,
                              int base_time_length, int fs, double window_length_in_time,
                              double *main_window) {
        double tmp = 0.0;
        for (int i = 0; i < base_time_length; ++i) {
            tmp = (index_raw[i] - 1.0) / fs - current_position;
            main_window[i] = 0.42 +
                             0.5 * cos(2.0 * Pi * tmp / window_length_in_time) +
                             0.08 * cos(4.0 * Pi * tmp / window_length_in_time);
        }
    }

//-----------------------------------------------------------------------------
// GetDiffWindow() generates the differentiated window.
// Diff means differential.
//-----------------------------------------------------------------------------
    static void GetDiffWindow(const double *main_window, int base_time_length,
                              double *diff_window) {
        diff_window[0] = -main_window[1] / 2.0;
        for (int i = 1; i < base_time_length - 1; ++i)
            diff_window[i] = -(main_window[i + 1] - main_window[i - 1]) / 2.0;
        diff_window[base_time_length - 1] = main_window[base_time_length - 2] / 2.0;
    }

//-----------------------------------------------------------------------------
// GetSpectra() calculates two spectra of the waveform windowed by windows
// (main window and diff window).
//-----------------------------------------------------------------------------
    static void GetSpectra(const double *x, int x_length, int fft_size,
                           const int *index_raw, const double *main_window, const double *diff_window,
                           int base_time_length, const ForwardRealFFT *forward_real_fft,
                           fft_complex *main_spectrum, fft_complex *diff_spectrum) {
        int *index = new int[base_time_length];

        for (int i = 0; i < base_time_length; ++i)
            index[i] = myMax<int>(0, myMin<int>(x_length - 1, index_raw[i] - 1));
        for (int i = 0; i < base_time_length; ++i)
            forward_real_fft->waveform[i] = x[index[i]] * main_window[i];
        for (int i = base_time_length; i < fft_size; ++i)
            forward_real_fft->waveform[i] = 0.0;

        fft_execute(forward_real_fft->forward_fft);
        for (int i = 0; i <= fft_size / 2; ++i) {
            main_spectrum[i][0] = forward_real_fft->spectrum[i][0];
            main_spectrum[i][1] = forward_real_fft->spectrum[i][1];
        }

        for (int i = 0; i < base_time_length; ++i)
            forward_real_fft->waveform[i] = x[index[i]] * diff_window[i];
        for (int i = base_time_length; i < fft_size; ++i)
            forward_real_fft->waveform[i] = 0.0;
        fft_execute(forward_real_fft->forward_fft);
        for (int i = 0; i <= fft_size / 2; ++i) {
            diff_spectrum[i][0] = forward_real_fft->spectrum[i][0];
            diff_spectrum[i][1] = forward_real_fft->spectrum[i][1];
        }

        delete[] index;
    }

//-----------------------------------------------------------------------------
// FixF0() fixed the F0 by instantaneous frequency.
//-----------------------------------------------------------------------------
    static double FixF0(const double *power_spectrum, const double *numerator_i,
                        int fft_size, int fs, double initial_f0, int number_of_harmonics) {
        double *amplitude_list = new double[number_of_harmonics];
        double *instantaneous_frequency_list = new double[number_of_harmonics];
        int index;
        for (int i = 0; i < number_of_harmonics; ++i) {
            index = matlab_round(initial_f0 * fft_size / fs * (i + 1));
            instantaneous_frequency_list[i] = power_spectrum[index] == 0.0 ? 0.0 :
                                              static_cast<double>(index) * fs / fft_size +
                                              numerator_i[index] / power_spectrum[index] * fs / 2.0 / Pi;
            amplitude_list[i] = sqrt(power_spectrum[index]);
        }
        double denominator = 0.0;
        double numerator = 0.0;
        for (int i = 0; i < number_of_harmonics; ++i) {
            numerator += amplitude_list[i] * instantaneous_frequency_list[i];
            denominator += amplitude_list[i] * (i + 1);
        }
        delete[] amplitude_list;
        delete[] instantaneous_frequency_list;
        return numerator / (denominator + SafeGuardMinimum);
    }

//-----------------------------------------------------------------------------
// GetTentativeF0() calculates the F0 based on the instantaneous frequency.
//-----------------------------------------------------------------------------
    static double GetTentativeF0(const double *power_spectrum,
                                 const double *numerator_i, int fft_size, int fs, double initial_f0) {
        double tentative_f0 = FixF0(power_spectrum, numerator_i, fft_size, fs, initial_f0, 2);

        // If the fixed value is too large, the result will be rejected.
        if (tentative_f0 <= 0.0 || tentative_f0 > initial_f0 * 2) return 0.0;

        return FixF0(power_spectrum, numerator_i, fft_size, fs, tentative_f0, 6);
    }

//-----------------------------------------------------------------------------
// GetMeanF0() calculates the instantaneous frequency.
//-----------------------------------------------------------------------------
    static double GetMeanF0(const double *x, int x_length, int fs,
                            double current_position, double initial_f0, int fft_size,
                            double window_length_in_time, const double *base_time,
                            int base_time_length) {
        ForwardRealFFT forward_real_fft = {0};
        InitializeForwardRealFFT(fft_size, &forward_real_fft);
        fft_complex *main_spectrum = new fft_complex[fft_size];
        fft_complex *diff_spectrum = new fft_complex[fft_size];

        int *index_raw = new int[base_time_length];
        double *main_window = new double[base_time_length];
        double *diff_window = new double[base_time_length];

        GetBaseIndex(current_position, base_time, base_time_length, fs, index_raw);
        GetMainWindow(current_position, index_raw, base_time_length, fs, window_length_in_time, main_window);
        GetDiffWindow(main_window, base_time_length, diff_window);
        GetSpectra(x, x_length, fft_size, index_raw, main_window, diff_window,
                   base_time_length, &forward_real_fft, main_spectrum, diff_spectrum);

        double *power_spectrum = new double[fft_size / 2 + 1];
        double *numerator_i = new double[fft_size / 2 + 1];
        for (int j = 0; j <= fft_size / 2; ++j) {
            numerator_i[j] = main_spectrum[j][0] * diff_spectrum[j][1] - main_spectrum[j][1] * diff_spectrum[j][0];
            power_spectrum[j] = main_spectrum[j][0] * main_spectrum[j][0] + main_spectrum[j][1] * main_spectrum[j][1];
        }

        double tentative_f0 = GetTentativeF0(power_spectrum, numerator_i, fft_size, fs, initial_f0);

        delete[] diff_spectrum;
        delete[] diff_window;
        delete[] main_window;
        delete[] index_raw;
        delete[] numerator_i;
        delete[] power_spectrum;
        delete[] main_spectrum;
        DestroyForwardRealFFT(&forward_real_fft);

        return tentative_f0;
    }

//-----------------------------------------------------------------------------
// GetRefinedF0() fixes the F0 estimated by Dio(). This function uses
// instantaneous frequency.
//-----------------------------------------------------------------------------
    static double GetRefinedF0(const double *x, int x_length, int fs,
                               double current_potision, double initial_f0) {
        if (initial_f0 <= FloorF0StoneMask || initial_f0 > fs / 12.0)
            return 0.0;

        int half_window_length = static_cast<int>(1.5 * fs / initial_f0 + 1.0);
        double window_length_in_time = (2.0 * half_window_length + 1.0) / fs;
        double *base_time = new double[half_window_length * 2 + 1];
        for (int i = 0; i < half_window_length * 2 + 1; i++)
            base_time[i] = static_cast<double>(-half_window_length + i) / fs;
        int fft_size = static_cast<int>(pow(2.0, 2.0 + (int)(log(half_window_length * 2.0 + 1.0) / Log2)));

        double mean_f0 = GetMeanF0(x, x_length, fs, current_potision,
                                   initial_f0, fft_size, window_length_in_time, base_time,
                                   half_window_length * 2 + 1);

        // If amount of correction is overlarge (20 %), initial F0 is employed.
        if (fabs(mean_f0 - initial_f0) / initial_f0 > 0.2) mean_f0 = initial_f0;

        delete[] base_time;

        return mean_f0;
    }

}  // namespace

EstimateF0WithDIO::EstimateF0WithDIO(double msFramePeriod)
        : msFramePeriod(msFramePeriod) {
}

int EstimateF0WithDIO::getF0LengthFor(unsigned int samplingFrequency, unsigned int waveLength) const {
    return GetSamplesForDIO(samplingFrequency, waveLength, msFramePeriod);
}

bool EstimateF0WithDIO::operator()(F0Contour *output, const Waveform *input) {
    if(output->length < getF0LengthFor(input->samplingFrequency, input->length)) {
        return false;
    }
    double *t = new double[output->length];
    const double channelsInOctave = 2.0;
    const int speed = 1;
    const double allowedRange = 0.1;
    DioGeneralBody(input->data, input->length, input->samplingFrequency, msFramePeriod,
            FloorF0, CeilF0, channelsInOctave, speed, allowedRange, t, output->data);
    for (int i = 0; i < output->length; i++)
        output->data[i] = GetRefinedF0(input->data, input->length, input->samplingFrequency, t[i], output->data[i]);
    delete[] t;
    return true;
}
