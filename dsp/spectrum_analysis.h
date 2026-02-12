#pragma once 

#include <math.h>

#include "../math/complex_math.h"

void calc_magnitudes(complex *dft_res, float *magnitudes, int N); 

void calc_phase_shifts(complex *dft_res, float *shifts, int N); 

void calc_power_spectrum(complex *dft_res, float *power_spectrum, int N); 

void calc_bin_frequencies(float *frequencies, int Fs, int N); 

void calc_bin_amplitudes(float *magnitudes, float *amplitudes, int N); 