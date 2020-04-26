// Modified by Hal@shurabaP.
// The original of this source code is in WORLD (https://github.com/mmorise/World/) by mmorise.
#ifndef __DSP_FFT_UTIL_HPP__
#define __DSP_FFT_UTIL_HPP__

// Complex number for FFT
typedef double fft_complex[2];
// Struct used for FFT
typedef struct {
    int n;
    int sign;
    unsigned int flags;
    fft_complex *c_in;
    double *in;
    fft_complex *c_out;
    double *out;
    double *input;
    int *ip;
    double *w;
} fft_plan;

// Commands for FFT (This is the same as FFTW)
#define FFT_FORWARD 1
#define FFT_BACKWARD 2
#define FFT_ESTIMATE 3

// Forward FFT in the real sequence
typedef struct {
    int fft_size;
    double *waveform;
    fft_complex *spectrum;
    fft_plan forward_fft;
} ForwardRealFFT;

// Inverse FFT in the real sequence
typedef struct {
    int fft_size;
    double *waveform;
    fft_complex *spectrum;
    fft_plan inverse_fft;
} InverseRealFFT;

// Inverse FFT in the complex sequence
typedef struct {
    int fft_size;
    fft_complex *input;
    fft_complex *output;
    fft_plan inverse_fft;
} InverseComplexFFT;

// Minimum phase analysis from logarithmic power spectrum
typedef struct {
    int fft_size;
    double *log_spectrum;
    fft_complex *minimum_phase_spectrum;
    fft_complex *cepstrum;
    fft_plan inverse_fft;
    fft_plan forward_fft;
} MinimumPhaseAnalysis;

//-----------------------------------------------------------------------------
// These functions are used to speed up the processing.
// Forward FFT
void InitializeForwardRealFFT(int fft_size, ForwardRealFFT *forward_real_fft);

void DestroyForwardRealFFT(ForwardRealFFT *forward_real_fft);

// Inverse FFT
void InitializeInverseRealFFT(int fft_size, InverseRealFFT *inverse_real_fft);

void DestroyInverseRealFFT(InverseRealFFT *inverse_real_fft);

// Inverse FFT (Complex)
void InitializeInverseComplexFFT(int fft_size,
                                 InverseComplexFFT *inverse_complex_fft);

void DestroyInverseComplexFFT(InverseComplexFFT *inverse_complex_fft);

// Minimum phase analysis (This analysis uses FFT)
void InitializeMinimumPhaseAnalysis(int fft_size,
                                    MinimumPhaseAnalysis *minimum_phase);

void GetMinimumPhaseSpectrum(const MinimumPhaseAnalysis *minimum_phase);

void DestroyMinimumPhaseAnalysis(MinimumPhaseAnalysis *minimum_phase);

#endif
