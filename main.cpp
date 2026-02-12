#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>


#include "math/complex_math.h"

#include "dsp/fourier_transform.h"
#include "dsp/spectrum_analysis.h"
#include "dsp/window_funcs.h"

#include "libs/matplotlibcpp.h"

#include "utils/memops.h"

namespace plt = matplotlibcpp;


typedef struct s_signal
{
	int N; 					// Number of DFT / IDFT points 

	int Fs; 				// Sample rate 

	WINDOW_T window; 	// Type of window to apply 
	
	float *input_signal;	// Pointer to input signal buffer 
	
	complex *dft_res;		// Pointer to DFT result 
	
	complex *idft_res;  	// Pointer to IDFT result 
	
	float *magnitudes;		// Pointer to calculated magnitudes
	
	float *phase_shifts;	// Pointer to calculated phase shifts 
	
	float *power_spectrum;	// Pointer to calculated power spectrum
	
	float *amplitudes; 		// Pointer to calculated ampllitudes
	
	float *frequencies;		// Pointer to calculated frequencies
	
} s_signal; 

void s_signal_init(s_signal *s, int N, int sample_rate, WINDOW_T win) 
{
	s->N = N; 

	s->Fs = sample_rate;  

	s->window = win; 
	
	s->input_signal = (float*)malloc_nc(N * sizeof(float)); 
	
	s->dft_res = (complex*)malloc_nc(N * sizeof(float)); 
	
	s->idft_res = (complex*)malloc_nc(N * sizeof(float)); 
	
	s->magnitudes = (float*)malloc_nc(N * sizeof(float));
	
	s->phase_shifts = (float*)malloc_nc(N * sizeof(float));
	
	s->power_spectrum = (float*)malloc_nc(N * sizeof(float));
	
	s->amplitudes = (float*)malloc_nc(N * sizeof(float));

	s->frequencies = (float*)malloc_nc(N * sizeof(float));
}

void s_signal_deinit(s_signal *s) 
{
	free_nc(s->input_signal); 
	
	free_nc(s->dft_res); 
	
	free_nc(s->idft_res); 
	
	free_nc(s->magnitudes);
	
	free_nc(s->phase_shifts);
	
	free_nc(s->power_spectrum);
	
	free_nc(s->amplitudes); 	
	
	free_nc(s->frequencies);
}

void display_peaks(float *x, float *y, int peaks_num) 
{
	std::vector<float> x_v(x, x + peaks_num);
	std::vector<float> y_v(y, y + peaks_num);

	// draw vertical line from 0 to y[i]
	for (size_t i = 0; i < x_v.size(); ++i) {
        plt::plot({x_v[i], x_v[i]}, {0, y_v[i]}, "b-"); // "b-" = blue line
    }

	// draw circles on peaks 
	plt::plot(x_v, y_v, "ro"); // "ro" = red circles 
    plt::show();
}

void analyze_audio(float *sound_data, s_signal *signal_str) 
{ 
	if(signal_str->window == NO_WINDOW)
	{ 
		// Calculate DFT and IDFT on non-windowed signal
		dft(sound_data, signal_str->dft_res, signal_str->N);
		idft(signal_str->dft_res, signal_str->idft_res, signal_str->N);
	}
	else 
	{ 
		// Apply signal Windowing to avoid spectral leakage 
		if(signal_str->window == HANN_WINDOW)
			hann_window(sound_data, signal_str->input_signal, signal_str->N);
	
		if(signal_str->window == HAMMING_WINDOW) 
			perror("Not implemented! \n"); 

		// Calculate DFT on a signal 
		dft(signal_str->input_signal, signal_str->dft_res, signal_str->N);
	}

	// Calculate magnitudes, phase shifts and power spectrum of a signal 
	calc_magnitudes(signal_str->dft_res, signal_str->magnitudes, signal_str->N); 
	calc_phase_shifts(signal_str->dft_res, signal_str->phase_shifts, signal_str->N); 
	calc_power_spectrum(signal_str->dft_res, signal_str->power_spectrum, signal_str->N); 

	float Fnyquist = (float)signal_str->Fs / 2; 			// Nyquist frequency 
	float delta_F = (float)signal_str->Fs / signal_str->N; 	// Frequency resolution - step between spectral bins  
	float T = (float)signal_str->N / signal_str->Fs; 		// Time window (time duration of analyzed signal)

	// Calculate frequency of each bin 
	calc_bin_frequencies(signal_str->frequencies, signal_str->Fs, signal_str->N); 

	// Calculate amplitude of each bin 
	calc_bin_amplitudes(signal_str->magnitudes, signal_str->amplitudes, signal_str->N); 
	
	int display_points = 256; 
	display_peaks(signal_str->frequencies, signal_str->amplitudes, display_points); 
}

void log_sf_info(SF_INFO *sf_info)
{ 
	printf("\n"); 
	printf("Sound file info: \n"); 
	printf("\t frames: %d\n", sf_info->frames); // total samples count in audiofile 
	printf("\t samplerate: %d\n", sf_info->samplerate); 
	printf("\t channels: %d\n", sf_info->channels); 
	printf("\t format: %d\n", sf_info->format); 
	printf("\t sections: %d\n", sf_info->sections); 
	printf("\t seekable: %d\n", sf_info->seekable); 

	float sf_time = (float)sf_info->frames / sf_info->samplerate; 
	printf("\t time: %.2f seconds\n", sf_time); 

	printf("\n\n"); 
}

SF_INFO sf_info;

s_signal signal_str; 

int main(int argc, char **argv)
{
	char sf_path[] = "data/ff-16b-2c-44100hz.wav"; 
	
	SNDFILE* f_sound = sf_open(sf_path, SFM_READ, &sf_info);

	log_sf_info(&sf_info); 

	int16_t *sound_data = (int16_t*)malloc_nc(sf_info.frames*sf_info.channels * sizeof(int16_t)); 
	int16_t *sound_data_ch0 = (int16_t*)malloc_nc(sf_info.frames*sf_info.channels/2 * sizeof(int16_t)); 
	int16_t *sound_data_ch1 = (int16_t*)malloc_nc(sf_info.frames/2 * sizeof(int16_t)); 
	float *sound_data_ch0_fl = (float*)malloc_nc(sf_info.frames*sf_info.channels/2 * sizeof(float)); 

	sf_read_short(f_sound, (int16_t*)sound_data, sf_info.frames*sf_info.channels);

	for(int i = 0, k = 0; i < sf_info.frames-2; k++, i += 2)
		sound_data_ch0[k] = sound_data[i];

	for(int i = 0; i < sf_info.frames/2; i++)
		sound_data_ch0_fl[i] = (float)sound_data_ch0[i] / INT16_MAX; 

	int N = 4096; 
	int sample_rate = sf_info.samplerate; 
	// s_signal_init(&signal_str, N, sample_rate, HANN_WINDOW); 
	s_signal_init(&signal_str, N, sample_rate, NO_WINDOW); 

	analyze_audio(sound_data_ch0_fl, &signal_str); 
	
	s_signal_deinit(&signal_str); 

	free_nc(sound_data); 

	return 0;  
}