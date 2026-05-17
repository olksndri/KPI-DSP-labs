#include <SDL2/SDL_pixels.h>
#include <cstdint>
#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>
#include <time.h>
#include <sstream>
#include <chrono>
#include <iostream>
#include <thread>
#include <iomanip>
#include <SDL2/SDL.h>
#include <SDL2/SDL_main.h>
#include "tiffio.h"
#include "image_compression/rle.h"
#include "image_compression/jpeg.h"

#include "math/complex_math.h"

#include "dsp/fourier_transform.h"
#include "dsp/spectrum_analysis.h"
#include "dsp/window_funcs.h"
#include "dsp/mel_filterbank.h"
#include "dsp/dct.h"

#include "libs/matplotlibcpp.h"

#include "utils/memops.h"

#include "utils/fileops.h"

#include "utils/timeops.h"

namespace plt = matplotlibcpp;


int use_logs = 0;
int use_time_profile = 0;

int lab = 4;


#include "raylib.h"
#include <string.h>
#include <stdint.h>

#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 600
#define HISTORY_SIZE SCREEN_WIDTH

float rawHistory[HISTORY_SIZE] = { 0 };
float procHistory[HISTORY_SIZE] = { 0 };


void tif_get_sizes(TIFF* f_tif, uint32_t &img_w, uint32_t &img_h, uint32_t &bitsize)
{
	TIFFGetMode(f_tif);
    TIFFGetField(f_tif, TIFFTAG_IMAGEWIDTH, &img_w);
    TIFFGetField(f_tif, TIFFTAG_IMAGELENGTH, &img_h);
    TIFFGetField(f_tif, TIFFTAG_BITSPERSAMPLE, &bitsize);
    printf("w %d\n", img_w);
	printf("h %d\n", img_h);
	printf("bitsize %d\n", bitsize);
}

void tif_read(TIFF* f_tif, uint32_t *image, uint32_t img_w, uint32_t img_h, uint32_t bitsize)
{
	if (f_tif == NULL)
	{
		printf("open image error\n");
		exit(0);
	}
	else
	{
	    if (image != NULL) {
	        // Read the image into the buffer
	        if (TIFFReadRGBAImageOriented(f_tif, img_w, img_h, image, 0)) {
	        } else {
	            // Handle error
				printf("read image error\n");
				exit(0);
	        }
	    }
	}
}

typedef struct s_signal
{
	int N; 					// Number of DFT / IDFT points

	int M; 					// Number of mel filterbank filters

	int Fs; 				// Sample rate

	WINDOW_T window; 	// Type of window to apply

	// float *input_signal;	// Pointer to input signal buffer
	complex *input_signal;	// Pointer to input signal buffer

	complex *dft_res;		// Pointer to DFT result

	complex *idft_res;  	// Pointer to IDFT result

	float *magnitudes;		// Pointer to calculated magnitudes

	float *phase_shifts;	// Pointer to calculated phase shifts

	float *power_spectrum;	// Pointer to calculated power spectrum

	float *amplitudes; 		// Pointer to calculated ampllitudes

	float *frequencies;		// Pointer to calculated frequencies

	float *mel_energies;	// Pointer to calculated mel energies

} s_signal;

void s_signal_init(s_signal *s, int N, int M, int sample_rate, WINDOW_T win)
{
	s->N = N;

	s->M = M;

	s->Fs = sample_rate;

	s->window = win;

	s->input_signal = (complex*)malloc_nc(N * sizeof(complex));

	s->dft_res = (complex*)malloc_nc(N * sizeof(complex));

	s->idft_res = (complex*)malloc_nc(N * sizeof(complex));

	s->magnitudes = (float*)malloc_nc(N * sizeof(float));

	s->phase_shifts = (float*)malloc_nc(N * sizeof(float));

	s->power_spectrum = (float*)malloc_nc(N * sizeof(float));

	s->amplitudes = (float*)malloc_nc(N * sizeof(float));

	s->frequencies = (float*)malloc_nc(N * sizeof(float));

	s->mel_energies = (float*)malloc_nc(M * sizeof(float));
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

	free_nc(s->mel_energies);
}

void display_peaks(float *x, float *y,
	const std::string plot_label, const std::string label_x, const std::string label_y, int peaks_num)
{

	// plt::figure_size(1200, 780);

	std::vector<float> x_v(x, x + peaks_num);
	std::vector<float> y_v(y, y + peaks_num);

	// draw vertical line from 0 to y[i]
	for (size_t i = 0; i < x_v.size(); ++i) {
        plt::plot({x_v[i], x_v[i]}, {0, y_v[i]}, "b-"); // "b-" = blue line
    }

	// draw circles on peaks
	plt::plot(x_v, y_v, "ro"); // "ro" = red circles

	plt::title(plot_label);
	plt::xlabel(label_x);
	plt::ylabel(label_y);

	// printf("showing graph\n");
    plt::show();
	// printf("ended showing graph\n");
    // plt::close();
}

void analyze_audio_0(float *sound_data, s_signal *signal_str)
{
	double t0, t1;

	float_to_complex(sound_data, signal_str->input_signal, signal_str->N);

	if(signal_str->window == NO_WINDOW)
	{
		// Calculate DFT and IDFT  on non-windowed signal
		t0 = current_time_ms();
		dft(signal_str->input_signal, signal_str->dft_res, signal_str->N);
		idft(signal_str->dft_res, signal_str->idft_res, signal_str->N);
		t1 = current_time_ms();
		printf("DFT + IDFT time: %.2lf [ms]\n", t1-t0);

		log_to_file("input_signal", signal_str->input_signal, signal_str->N, COMPLEX);
		log_to_file("dft_res", signal_str->dft_res, signal_str->N, COMPLEX);
		log_to_file("idft_res", signal_str->idft_res, signal_str->N, COMPLEX);

		t0 = current_time_ms();
		fft_recursive(signal_str->input_signal, signal_str->dft_res, signal_str->N);
		ifft_recursive(signal_str->dft_res, signal_str->idft_res, signal_str->N);
		ifft_recursive_scale(signal_str->idft_res, signal_str->N);
		t1 = current_time_ms();
		printf("DFFT + IDFFT time: %.2lf [ms]\n", t1-t0);

		log_to_file("dfft_recursive_res", signal_str->dft_res, signal_str->N, COMPLEX);
		log_to_file("idfft_recursive_res", signal_str->idft_res, signal_str->N, COMPLEX);
	}
	else
	{
		// Apply signal Windowing to avoid spectral leakage
		if(signal_str->window == HANN_WINDOW)
			hann_window_complex(signal_str->input_signal, signal_str->input_signal, signal_str->N);
		else
			perror("Not implemented! \n");

		// Calculate DFT on a signal
		t0 = current_time_ms();
		dft(signal_str->input_signal, signal_str->dft_res, signal_str->N);
		t1 = current_time_ms();
		printf("DFT time: %.2lf [ms]\n", t1-t0);

		t0 = current_time_ms();
		fft_recursive(signal_str->input_signal, signal_str->dft_res, signal_str->N);
		t1 = current_time_ms();
		printf("FFT time: %.2lf [ms]\n", t1-t0);
	}

	// Calculate magnitudes, phase shifts and power spectrum of a signal
	calc_magnitudes(signal_str->dft_res, signal_str->magnitudes, signal_str->N);
	calc_phase_shifts(signal_str->dft_res, signal_str->phase_shifts, signal_str->N);
	calc_power_spectrum(signal_str->dft_res, signal_str->power_spectrum, signal_str->N);

	float Fnyquist = (float)signal_str->Fs / 2; 			// Nyquist frequency [Hz]
	float delta_F = (float)signal_str->Fs / signal_str->N; 	// Frequency resolution - step between spectral bins, [Hz]
	float T = (float)signal_str->N / signal_str->Fs; 		// Time window (time duration of analyzed signal interval), [s]
	float T_ms = T * 1000 ; // Time window in [ms]

	printf("Nyquist frequency: %.2f [Hz]\n", Fnyquist);
	printf("Frequency resolution: %.2f [Hz]\n", delta_F);
	printf("Analysis time window: %.2f [ms]\n", T_ms);

	// Calculate frequency of each bin
	calc_bin_frequencies(signal_str->frequencies, signal_str->Fs, signal_str->N);

	// Calculate amplitude of each bin
	calc_bin_amplitudes(signal_str->magnitudes, signal_str->amplitudes, signal_str->N);

	int display_points = signal_str->N;
	std::string plot_title = "Plot title" ;
	display_peaks(signal_str->frequencies, signal_str->amplitudes, plot_title, "Frequency", "Amplitude", display_points);
}


// LAB 1
// Розробка програмного забезпечення для швидкого
// перетворення Фур’є голосового сигналу.

// Завдання: розробити програмне забезпечення для
// зчитування аудіоданних із звукового
// файлу, їх фрагментації та швидкого перетворення Фур’є кожного із
// фрагментів.
void analyze_audio_1(float *sound_data, s_signal *signal_str, int use_fft, int frame_num)
{
	double t0, t1;

	float_to_complex(sound_data, signal_str->input_signal, signal_str->N);

	// Apply signal Windowing to avoid spectral leakage
	if(signal_str->window == HANN_WINDOW)
		hann_window_complex(signal_str->input_signal, signal_str->input_signal, signal_str->N);

	// Calculate DFT and IDFT  on non-windowed signal
	t0 = current_time_ms();
	if(use_fft)
		fft_recursive(signal_str->input_signal, signal_str->dft_res, signal_str->N);
	else
		dft(signal_str->input_signal, signal_str->dft_res, signal_str->N);
	t1 = current_time_ms();

	printf("Fourier transform time: %.2lf [ms]\n", t1-t0);

	// Calculate magnitudes of a signal
	calc_magnitudes(signal_str->dft_res, signal_str->magnitudes, signal_str->N);

	float Fnyquist = (float)signal_str->Fs / 2; 			// Nyquist frequency [Hz]
	float delta_F = (float)signal_str->Fs / signal_str->N; 	// Frequency resolution - step between spectral bins, [Hz]
	float T = (float)signal_str->N / signal_str->Fs; 		// Time window (time duration of analyzed signal interval), [s]
	float T_ms = T * 1000 ; // Time window in [ms]

	printf("Nyquist frequency: %.2f [Hz]\n", Fnyquist);
	printf("Frequency resolution: %.2f [Hz]\n", delta_F);
	printf("Analysis time window: %.2f [ms]\n", T_ms);

	// Calculate frequency and amplitude of each bin
	calc_bin_frequencies(signal_str->frequencies, signal_str->Fs, signal_str->N);
	calc_bin_amplitudes(signal_str->magnitudes, signal_str->amplitudes, signal_str->N);


	std::stringstream ss;

	float time_window_beg = frame_num * T_ms;
	float time_window_end = time_window_beg + T_ms;
	std::string time_mes_unit = "[ms]";

	if(time_window_beg >= 1000)
	{
		time_window_beg /= 1000;
		time_window_end /= 1000;
		time_mes_unit = "[s]";
	}
	ss << std::fixed << std::setprecision(2);
	ss << "Analyzed time window: " << time_window_beg << "..." << time_window_end << " " << time_mes_unit;

	std::string plot_title = ss.str();

	int display_points = signal_str->N;
	display_peaks(signal_str->frequencies, signal_str->amplitudes, plot_title, "Frequency", "Amplitude", display_points);
}

FILE *gp;

void plot_current_mfcc(float* mfcc, int M, int frame_idx)
{
    fprintf(gp, "set title 'MFCC Frame: %d'\n", frame_idx);
    fprintf(gp, "set ylabel 'Value (dB/Amplitude)'\n");
    fprintf(gp, "set xlabel 'Coefficient Index'\n");

    fprintf(gp, "set style fill solid 1.0 border -1\n");
    fprintf(gp, "set boxwidth 0.8\n");
    fprintf(gp, "set grid y\n");

    fprintf(gp, "set autoscale y\n");
    fprintf(gp, "set autoscale x\n");

    fprintf(gp, "plot '-' with boxes title 'MFCC Coefficients'\n");

    for (int m = 0; m < M; m++) {
        fprintf(gp, "%d %f\n", m, mfcc[m]);
    }

    fprintf(gp, "e\n");
    fflush(gp);

    // pclose(gp);
}

void wait_for_1_key() {
    int x = 0;
    while (x != 1)
        std::cin >> x;
}

// LAB 2
// Розробка програмного забезпечення для розрахунку
// мел-кепстральних коефіцієнтів голосового сигналу.

// Завдання:
// розробити програмне забезпечення для розрахунку мел-кепстральних коефіцієнтів
// стаціонарних фрагментів голосового сигналу.
void calc_mfcc(float *sound_data, float *mel_filterbank, float *scratch, s_signal *signal_str, int samples)
{
    // Next steps need to be followed to calculate MFCCs:

    // 1. Pre-emphasize the signal: Amplify higher frequencies to balance the spectrum.
    float alpha = 0.96; // filter coefficient
    pre_emphasis_filter(scratch, sound_data, samples, alpha);

    float *mfcc = (float *)malloc_nc(signal_str->M * sizeof(float));

    int frame_count = 0;

    // 2. Framing: Break the signal into small, overlapping frames.
    int stride = signal_str->N / 2; // 50% overlap
    for(int i = 0; i <= samples - signal_str->N; i += stride)
    {
        float *cur_frame = scratch + i;

        float_to_complex(cur_frame, signal_str->input_signal, signal_str->N);

        // 3. Windowing: To soften the edges of each frame, apply a Hann window.
        hann_window_complex(signal_str->input_signal, signal_str->input_signal, signal_str->N);

       	// 4. FFT: Convert each frame from the time domain to the frequency domain.
        fft_recursive(signal_str->input_signal, signal_str->dft_res, signal_str->N);

        // 5. Mel-filterbank: Apply overlapping triangular filters spaced according to the Mel-scale on signals power spectrum.
        calc_power_spectrum(signal_str->dft_res, signal_str->power_spectrum, signal_str->N);

        apply_mel_filterbank(signal_str->mel_energies, signal_str->power_spectrum, mel_filterbank, signal_str->M, signal_str->N);

        // 6. Logarithm: To replicate the way a human ear reacts to sound strength take the logarithm of the filterbank outputs.
        for(int m = 0; m < signal_str->M; m++) // Convert mel energies into dB
            signal_str->mel_energies[m] = 10 * log10f(signal_str->mel_energies[m] + 1e-10);

        // 7. DCT: Apply the DCT to the log Mel-spectrum to obtain the Mel-frequency Cepstral Coefficients.
        dct_1d(mfcc, signal_str->mel_energies, signal_str->M);

        // 8. Visualize current frame MFCC using GNUplot
        gp = popen("gnuplot", "w");
        if (!gp) return;
        plot_current_mfcc(mfcc, signal_str->M, frame_count);

        // 9. Pause until user enters "1" in console
        wait_for_1_key();
        pclose(gp);
        frame_count++;
    }

    free_nc(mfcc);
}


// Exponential moving average filter
float ema(float yn_prev, float xn, float alpha)
{
	float yn = (1 - alpha) * yn_prev + alpha * xn;

	return yn;
}

// // LAB 3
// // Розробка програмного забезпечення для визначення меж
// // окремих слів в голосовому сигналі.

// // Завдання:
// // розробити програмне забезпечення для визначення часових кордонів пауз та
void vad(float *sound_data, float *mel_filterbank, float *scratch, s_signal *signal_str,
         int samples, Wave &wave, Sound &sound, float *rawHistory, float *procHistory)
{
    float pre_emphasis_alpha = 0.96;
    pre_emphasis_filter(scratch, sound_data, samples, pre_emphasis_alpha);

    float *mfcc = (float *)malloc_nc(signal_str->M * sizeof(float));
    int frame_count = 0;
    float energy_prev = 0.0f;
    float energy_noise = 0.0f;

    const int WINDOW_SIZE = signal_str->N;
    const int INFERENCES_SIZE = 5;

    int inferences[INFERENCES_SIZE];
    memset(inferences, 0, INFERENCES_SIZE * sizeof(int));

    int infer_cnt = 0;

    int is_voice = 0;
    int currentSample = 0;
    float elapsedSeconds = 0.0f;

    const int thresh_dB = 10.0f;

    int not_voice_cnt = 0;

    memset(rawHistory, 0, sizeof(float) * SCREEN_WIDTH);
    memset(procHistory, 0, sizeof(float) * SCREEN_WIDTH);

    while (!WindowShouldClose()) {
        if (IsSoundPlaying(sound)) {
            elapsedSeconds += GetFrameTime();
            long targetSample = (long)(elapsedSeconds * wave.sampleRate);

            while (currentSample + WINDOW_SIZE < targetSample &&
                   currentSample + WINDOW_SIZE < samples)
            {
            	float* current_window = &sound_data[currentSample];
                float* current_window_filtered = &scratch[currentSample];

                // --- MFFC --- //
                float_to_complex(current_window_filtered, signal_str->input_signal, signal_str->N);
                hann_window_complex(signal_str->input_signal, signal_str->input_signal, signal_str->N);
                fft_recursive(signal_str->input_signal, signal_str->dft_res, signal_str->N);
                calc_power_spectrum(signal_str->dft_res, signal_str->power_spectrum, signal_str->N);
                apply_mel_filterbank(signal_str->mel_energies, signal_str->power_spectrum, mel_filterbank, signal_str->M, signal_str->N);

                for(int m = 0; m < signal_str->M; m++)
                    signal_str->mel_energies[m] = 10 * log10f(signal_str->mel_energies[m] + 1e-10);

                dct_1d(mfcc, signal_str->mel_energies, signal_str->M);

                // --- VAD logic --- //
                float energy = mfcc[0];

                // Exponential moving average smoothing
                energy = ema(energy_prev, energy, 0.45);
                is_voice = (energy > (energy_noise + thresh_dB)) ? 1 : 0;
                energy_prev = energy;

                inferences[infer_cnt] = is_voice;
                infer_cnt++;

                // Dynamic noise floor tracking
                int is_noise_sequence = 1;
                for(int i = 0; i < INFERENCES_SIZE; i++)
                {
                	if(inferences[i])
                	{
                 		is_noise_sequence = 0;
                  		break;
                 	}
                }

                if(is_noise_sequence)
                	energy_noise = energy;

                if(infer_cnt == INFERENCES_SIZE)
               		infer_cnt = 0;

                // --- Visualization ---
                memmove(rawHistory, &rawHistory[1], (SCREEN_WIDTH - 1) * sizeof(float));
                memmove(procHistory, &procHistory[1], (SCREEN_WIDTH - 1) * sizeof(float));

                rawHistory[SCREEN_WIDTH - 1] = current_window[0];
                procHistory[SCREEN_WIDTH - 1] = is_voice ? 0.8f : 0.0f;

                currentSample += WINDOW_SIZE;

                frame_count++;
            }
        }

        // --- Drawing ---
        BeginDrawing();
            ClearBackground(BLACK);

            DrawText("Original Signal", 10, 10, 20, RAYWHITE);
            for (int x = 0; x < SCREEN_WIDTH; x++) {
                DrawLine(x, 150, x, 150 + (int)(rawHistory[x] * 200), BLUE);
            }

            DrawLine(0, 300, SCREEN_WIDTH, 300, DARKGRAY);

            DrawText("VAD Decision (Green = Speech)", 10, 310, 20, RAYWHITE);
            for (int x = 0; x < SCREEN_WIDTH; x++) {
                if (procHistory[x] > 0.1f) {
                    DrawLine(x, 350, x, 550, GREEN);
                }
            }
            DrawFPS(10, 10);
        EndDrawing();
    }
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

s_signal signal;

int main(int argc, char **argv)
{
	try {
        plt::backend("Qt5Agg");
    } catch (const std::exception& e) {
        printf("Backend error: %s\n", e.what());
    }


	if(lab == 1 || lab == 2)
	{
		SNDFILE* f_sound;
		if(argc > 1)
		{
 		f_sound = sf_open(argv[1], SFM_READ, &sf_info);

			if(f_sound == NULL)
			{
				perror("Error occured while opening file: ");
				printf("%s\n", argv[1]);
				abort();
			}
		}
		else
		{
			perror("Please provide path to sound file to be analyzed!\n");
			abort();
		}

		log_sf_info(&sf_info);
		int total_samples = sf_info.frames * sf_info.channels;

		int16_t *sound_data = (int16_t*)malloc_nc(total_samples * sizeof(int16_t));
		int16_t *sound_data_ch0 = (int16_t*)malloc_nc(sf_info.frames * sizeof(int16_t));
		float *sound_data_ch0_fl = (float*)malloc_nc(sf_info.frames * sizeof(float));
		float *scratch = (float*)malloc_nc(sf_info.frames * sizeof(float));

		sf_read_short(f_sound, (int16_t*)sound_data, sf_info.frames*sf_info.channels);

		for(int i = 0, k = 0; i < total_samples; k++, i += sf_info.channels)
			sound_data_ch0[k] = sound_data[i];

		for(int i = 0; i < sf_info.frames; i++)
			sound_data_ch0_fl[i] = (float)sound_data_ch0[i] / INT16_MAX;


		int sample_rate = sf_info.samplerate;

     	if (lab == 1)
      	{
      		int N = 0;
       		int M = 0;

         	int win_flag = 0;
         	int use_fft = 0;

         	std::cout << "Enter number of DFT points: \t";
         	std::cin >> N;
         	std::cout << "\n";
         	std::cout << "Should we use FFT? Print 1 or 0: \t";
         	std::cin >> use_fft;
         	std::cout << "\n";
         	std::cout << "Should we use hann window? Print 1 or 0: \t";
         	std::cin >> win_flag;
         	std::cout << "\n";

         	printf("Number of DFT points: %d\n", N);
         	printf("Used window: %s\n", (win_flag) ? "HANN_WINDOW" :"NO_WINDOW");

         	s_signal_init(&signal, N, M, sample_rate, (win_flag) ? HANN_WINDOW : NO_WINDOW);
         	int frames_to_process = sf_info.frames / N;
         	for(int i = 0; i < frames_to_process; i++)
         	{
         		analyze_audio_1(&sound_data_ch0_fl[N*i], &signal, use_fft, i);
         	}

         	s_signal_deinit(&signal);
       	}
		else
		{
			int N = 1024; // FFT pointsя
			int M = 40; // filters number
			int Fs = sample_rate;

			s_signal_init(&signal, N, M, sample_rate, HANN_WINDOW);

			float* mel_filterbank = create_mel_filterbank(Fs, N, M);

			calc_mfcc(sound_data_ch0_fl, mel_filterbank, scratch, &signal, sf_info.frames);

			destroy_mel_filterbank(mel_filterbank);

			s_signal_deinit(&signal);
		}

     	free_nc(scratch);
		free_nc(sound_data_ch0_fl);
		free_nc(sound_data_ch0);
		free_nc(sound_data);

		sf_close(f_sound);
	}
	else if (lab == 3)
	{
	 	InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Real-time Audio Scroll");
		InitAudioDevice();

		Wave wave = LoadWave(argv[1]);
		Sound sound = LoadSoundFromWave(wave);

		int currentSample = 0;

		PlaySound(sound);
		SetTargetFPS(25);

		WaveFormat(&wave, wave.sampleRate, 32, 1);

		float *samples = (float*)malloc_nc(wave.frameCount * sizeof(float));
		memcpy(samples, wave.data, wave.frameCount * sizeof(float));

		float *scratch = (float*)malloc_nc(wave.frameCount * sizeof(float));

		int N = 1024; // FFT points
		int M = 40; // filters number
		int Fs = wave.sampleRate;

		s_signal_init(&signal, N, M, Fs, HANN_WINDOW);

		float* mel_filterbank = create_mel_filterbank(Fs, N, M);

		vad(samples, mel_filterbank, scratch, &signal, wave.frameCount / wave.channels, wave, sound, rawHistory, procHistory);

		destroy_mel_filterbank(mel_filterbank);

		s_signal_deinit(&signal);

		free_nc(scratch);
	}
	else if (lab == 4)
	{
		TIFF* f_tif = TIFFOpen(argv[1], "r");

		uint32_t bitsize;
		uint32_t img_w, img_h;
		tif_get_sizes(f_tif, img_w, img_h, bitsize);

		uint32_t *image = (uint32_t*) _TIFFmalloc(img_h * img_w * sizeof(uint32_t));
		tif_read(f_tif, image, img_w, img_h, bitsize);

		const int y_margin = 100;
		const int x_margin = 25;

		const int window_w = (img_w * 2) + x_margin * 4 + x_margin;
		const int window_h = (img_h) + y_margin * 2;

		SDL_Window* window = NULL;
	    SDL_Renderer* renderer = NULL;
		SDL_Texture* texture = NULL;

	    // Initialize SDL
	    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
	        fprintf(stderr, "SDL could not be initialized! SDL_Error: %s\n", SDL_GetError());
	        return 1;
	    }

	    // Create window and renderer
	    if (SDL_CreateWindowAndRenderer(window_w, window_h, SDL_WINDOW_SHOWN, &window, &renderer) < 0) {
	        fprintf(stderr, "Window and renderer could not be created! SDL_Error: %s\n", SDL_GetError());
	        return 1;
	    }

		texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, window_w, window_h);

		uint32_t window_buff[window_w * window_h];

		// SDL_Color black = (SDL_Color) {.r = 0, .g = 0, .b = 0, .a = 255};
		// SDL_Color *p_win = (SDL_Color *)window_buff;
		// for(int i = 0; i < window_h * window_w; i++ )
		// 	p_win[i] = black;

		for(int y = 0; y < img_h; y++)
		{
			for(int x = 0; x < img_w; x++)
			{
				int inp_idx = (y * img_w) + x;
				int out_idx = ((y + y_margin) * window_w) + x + x_margin;
				window_buff[out_idx] = image[inp_idx];
			}
		}

		int orig_img_size = img_h * img_w * 4;
		// RLE_IMAGE rle_image_st;
		uint32_t* decompressed_image = (uint32_t* )malloc_nc(img_h * img_w * sizeof(uint32_t));
		uint8_t* RGB_in = (uint8_t* )malloc_nc(img_h * img_w * 3 * sizeof(uint8_t));
		uint8_t* RGB_out = (uint8_t* )malloc_nc(img_h * img_w * 3 * sizeof(uint8_t));
		// rle_compress(image, rle_image_st, img_w, img_h);
		// rle_decompress(decompressed_image, rle_image_st);

		for(int i = 0; i < img_h * img_w; i++)
		{
			uint32_t v = image[i];
			uint8_t r = v >> 16;
			uint8_t g = v >> 8;
			uint8_t b = v >> 0;
			RGB_in[i*3+0] = r;
			RGB_in[i*3+1] = g;
			RGB_in[i*3+2] = b;
		}

		PADDING pad_luminance;
		PADDING pad_chroma;

		RLE_IMAGE_1CH_EOB Y_rle;
		RLE_IMAGE_1CH_EOB Cb_rle;
		RLE_IMAGE_1CH_EOB Cr_rle;
		int use_chroma_downsampling = 1;
		jpeg_encode(Y_rle, Cb_rle, Cr_rle, pad_luminance, pad_chroma, RGB_in, img_h, img_w, use_chroma_downsampling);
		jpeg_decode(RGB_out, Y_rle, Cb_rle, Cr_rle, pad_luminance, pad_chroma, img_h, img_w, use_chroma_downsampling);

		for(int i = 0; i < img_h * img_w; i++)
		{
			uint8_t r = RGB_out[i*3+0];
			uint8_t g = RGB_out[i*3+1];
			uint8_t b = RGB_out[i*3+2];
			uint32_t rgba = (255 << 24) | (r << 16) | (g << 8) | (b << 0);
			decompressed_image[i]  = rgba ;
		}

		int compressed_img_data_size_bytes
			=	(Y_rle.data.size() * sizeof(int16_t)) +
				(Cb_rle.data.size() * sizeof(int16_t)) +
				(Cr_rle.data.size() * sizeof(int16_t));
		int compressed_img_metadata_size_bytes
			= 	(Y_rle.EOB.size() + sizeof(uint8_t)) +
				(Cb_rle.EOB.size() + sizeof(uint8_t)) +
				(Cr_rle.EOB.size() + sizeof(uint8_t)) +
				sizeof(pad_luminance) + sizeof(pad_chroma) +
				sizeof(img_h) + sizeof(img_w);

		int compressed_img_size = compressed_img_data_size_bytes + compressed_img_metadata_size_bytes;

		printf("original image size: %d bytes\n", orig_img_size);
		printf("compessed image size: %d bytes\n", compressed_img_size);

		for(int y = 0; y < img_h; y++)
		{
			for(int x = 0; x < img_w; x++)
			{
				int inp_idx = (y * img_w) + x;
				int out_idx = ((y + y_margin) * window_w) + x + 2 * x_margin + img_w;
				window_buff[out_idx] = decompressed_image[inp_idx];
			}
		}

	    // Main loop flag
	    int quit = 0;
	    // Event handler
	    SDL_Event e;

	    // Main loop
	    while (!quit) {
	        // Handle events on queue
	        while (SDL_PollEvent(&e) != 0) {
	            // User requests quit
	            if (e.type == SDL_QUIT) {
	                quit = 1;
	            }
	        }

	        // Clear screen with a color (e.g., black)
			// RGBA format
			SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
	        SDL_RenderClear(renderer);
			// memcpy(renderer, f_tif, sizeof(uint8_t) * 512 ) ;
			SDL_UpdateTexture(texture, NULL, window_buff, window_w * sizeof(uint8_t) * 4);
   			SDL_RenderCopy(renderer, texture, NULL, NULL); // Copy the texture to the renderer

      		// Update the screen
	        SDL_RenderPresent(renderer);
	    }

	    // Clean up
	    SDL_DestroyRenderer(renderer);
	    SDL_DestroyWindow(window);
	    SDL_Quit();

		_TIFFfree(image);
		TIFFClose(f_tif);
	}
	else
	{
	    std::cout   << "Error, there is no lab " << lab << "implementation."
					<< "Exiting program\n\n";
		exit(-1);
	}

	return 0;
}
