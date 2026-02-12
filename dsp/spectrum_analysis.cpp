#include "spectrum_analysis.h"

void calc_magnitudes(complex *dft_res, float *magnitudes, int N)
{
	for(int m = 0; m < N; m++)
	{ 
		complex v = dft_res[m];
		float mag = sqrtf(v.real*v.real + v.imag*v.imag); 
		magnitudes[m] = mag; 
	}
}

void calc_phase_shifts(complex *dft_res, float *shifts, int N)
{
	for(int m = 0; m < N; m++)
	{ 
		complex v = dft_res[m];
		float shift = atan2f(v.imag, v.real);
		shifts[m] = shift; 
	}
}

void calc_power_spectrum(complex *dft_res, float *power_spectrum, int N) 
{ 
	for(int m = 0; m < N; m++)
	{ 
		complex v = dft_res[m];
		float power = v.real*v.real + v.imag*v.imag; 
		power_spectrum[m] = power; 
	}
}

void calc_bin_frequencies(float *frequencies, int Fs, int N)
{ 	
	// Real signal has unique frequencies up to N / 2 
	// Above the N / 2 frequencies are "mirrored"
	int kmax = N / 2; 

	// Frequency resolution - step between spectral bins  
	float delta_F = (float)Fs / N; 

	for(int k = 0; k <= kmax; k++)
		frequencies[k] = k * delta_F; 
}


void calc_bin_amplitudes(float *magnitudes, float *amplitudes,  int N)
{ 
	// Real signal has unique frequencies up to N / 2 
	// Above N / 2 frequencies are "mirrored"
	int kmax = N / 2; 

	// DC amplitude 
	amplitudes[0] = magnitudes[0] / N; 

	for(int k = 1; k < kmax; k++)
		amplitudes[k] = (2.0f * magnitudes[k]) / N; 

	// Nyquist frequency amplitude 
	amplitudes[kmax] = magnitudes[kmax] / N; 
}
