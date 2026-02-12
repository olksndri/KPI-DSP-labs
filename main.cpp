#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>

#if 1 
#include <SDL2/SDL.h> // for initializing and shutdown functions
#include <SDL2/SDL_image.h> // for rendering images and graphics on screen
#include <SDL2/SDL_timer.h> // for using SDL_Delay() functions
#endif 

#include "include/t.h"


#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

// #include "gnuplot_i.h" // Include the library header



// Розробка програмного забезпечення для швидкого
// перетворення Фур’є голосового сигналу.

// Завдання: розробити програмне забезпечення для
// зчитування аудіоданних із звукового 
// файлу, їх фрагментації та швидкого перетворення Фур’є кожного із
// фрагментів.

typedef struct complex {
	float real; 
	float imag; 
} complex; 

complex complex_sub(complex a, complex b) 
{
	complex res = {
		.real = a.real - b.real, 
		.imag = a.imag - b.imag,
	}; 

	return res; 
}

complex complex_add(complex a, complex b) 
{
	complex res = {
		.real = a.real + b.real, 
		.imag = a.imag + b.imag,
	}; 

	return res; 
}

complex complex_mul(complex a, complex b) 
{
	complex res = {
		.real = a.real*b.real - a.imag*b.imag, 
		.imag = a.real*b.imag + b.real*a.imag,
	}; 

	return res; 
}

void dft(float *in, complex *out, int N) 
{
	for(int k = 0; k < N; k++)
	{ 
		// float re = 0; 
		// float im = 0; 
		complex Xout = { 0, 0 }; 

		for(int n = 0; n < N; n++)
		{
			float Xin = in[n];

			Xout.real += cosf((2*M_PI*n*k)/N) * Xin; 
			Xout.imag += -sinf((2*M_PI*n*k)/N) * Xin; 
		}
		
		// complex Xout = { .real = re, .imag = im }; 
		
		out[k] = Xout; 
	}
}


void idft(complex *in, complex *out, int N) 
{
	for(int k = 0; k < N; k++)
	{ 
		complex Xout = { .real = 0, .imag = 0 }; 

		for(int n = 0; n < N; n++)
		{
			complex Xin = in[n];
			
			complex basis_func = { 
				.real = cosf((2*M_PI*n*k)/N),
				.imag = sinf((2*M_PI*n*k)/N), 
			}; 

			Xout = complex_add(Xout, complex_mul(Xin, basis_func));  
		}

		Xout.real /= N; 
		Xout.imag /= N; 

		out[k] = Xout; 
	}
}

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
		// float im = dft_res[m].imag; 
		// float mag = magnitudes[m]; 
		// float shift = atanf(im/mag); 
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

void calc_bins_frequencies(float *frequencies, int Fs, int N)
{ 	
	// Real signal has unique frequencies up to N / 2 
	// Above the N / 2 frequencies are "mirrored"
	int kmax = N / 2; 

	// Frequency resolution - step between spectral bins  
	float delta_F = (float)Fs / N; 

	for(int k = 0; k <= kmax; k++)
		frequencies[k] = k * delta_F; 
}


void calc_bins_amplitudes(float *magnitudes, float *amplitudes,  int N)
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

void hann_window(float *data, float *windowed_data, int N)
{ 

	for(int n = 0; n < N; n++)
	{
		float window = 0.5 * (1 - cosf((2*M_PI*n)/(N-1))); 
		windowed_data[n] = data[n] * window; 
	}
}


typedef enum window_type 
{ 
	NO_WINDOW,
	HANN_WINDOW, 
	HAMMING_WINDOW, 
} window_type; 


void* malloc_nc(size_t size) 
{
    void *ptr = malloc(size);  

	if(ptr == NULL)
	{ 
		printf("Error in memory allocation! Aborting program... \n"); 
		abort();
	} 

	return ptr; 
}

void free_nc(void *ptr)
{
	if(ptr != NULL)
		free(ptr); 
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

int main()
{
	char sf_path[] = "data/ff-16b-2c-44100hz.wav"; 
	
	SNDFILE* f_sound = sf_open(sf_path, SFM_READ, &sf_info);

	log_sf_info(&sf_info); 

	// float *sound_data = (float*)malloc_nc(sf_info.frames * sizeof(float)); 
	int16_t *sound_data = (int16_t*)malloc_nc(sf_info.frames*sf_info.channels * sizeof(int16_t)); 
	int16_t *sound_data_ch0 = (int16_t*)malloc_nc(sf_info.frames*sf_info.channels/2 * sizeof(int16_t)); 
	float *sound_data_ch0_fl = (float*)malloc_nc(sf_info.frames*sf_info.channels/2 * sizeof(float)); 
	// uint16_t *sound_data_ch1 = (uint16_t*)malloc_nc(sf_info.frames/2 * sizeof(uint16_t)); 

	// sf_count_t read_data = 
	// sf_read_float(f_sound, sound_data, sf_info.frames) ;
	sf_read_short(f_sound, (short*)sound_data, sf_info.frames) ;
	// printf("read data: %d\n", read_data); 

	for(int i = 0, k = 0; i < sf_info.frames-2; k++, i += 2)
		sound_data_ch0[k] = sound_data[i];

	for(int i = 0; i < sf_info.frames/2; i++)
		sound_data_ch0_fl[i] = (float)sound_data_ch0[i] / INT16_MAX; 

	int sample_rate = sf_info.samplerate; 

// DFT/IDFT size 
// #define N 2048
	#define N 4096
	float windowed_sound[N]; 
	complex dft_arr[N]; 
	complex dft_windowed_arr[N]; 
	complex idft_arr[N];  
	float magnitudes[N];
	float shifts[N];
	float power_spectrum[N];
	float frequencies[N];
	float amplitudes[N];

	int ptr_shift = N; 

	hann_window(sound_data_ch0_fl+ptr_shift, windowed_sound, N);
	dft(sound_data_ch0_fl+ptr_shift, dft_arr, N);
	idft(dft_arr, idft_arr, N);

	dft(windowed_sound, dft_windowed_arr, N);
	calc_magnitudes(dft_windowed_arr, magnitudes, N); 
	calc_phase_shifts(dft_windowed_arr, shifts, N); 
	calc_power_spectrum(dft_windowed_arr, power_spectrum, N); 

	// Sample Rate 
	int Fs = sf_info.samplerate; 

	// Nyquist frequency 
	float Fnyquist = (float)Fs / 2; 

	// Frequency resolution - step between spectral bins  
	float delta_F = (float)Fs / N; 

	// Time window (length of analyzed signal)
	float T = (float)N / Fs; 

	calc_bins_frequencies(frequencies, Fs, N); 

	calc_bins_amplitudes(magnitudes, amplitudes, N); 

	
	// Frequencies above N / 2 are mirrored in real signal 
	int half_N = N / 2 + 1;
	
	int display_scale = 128; 

	
	// std::vector<float> freq_vec(frequencies, frequencies+half_N);
	// std::vector<float> ampl_vec(amplitudes, amplitudes+half_N);

	std::vector<float> freq_vec(frequencies, frequencies+display_scale);
	std::vector<float> ampl_vec(amplitudes, amplitudes+display_scale);


    for (size_t i = 0; i < freq_vec.size(); ++i) {
        // малюємо вертикальну лінію від 0 до y[i]
        plt::plot({freq_vec[i], freq_vec[i]}, {0, ampl_vec[i]}, "b-"); // "b-" = синя лінія
    }
	// можна додати кружечки на піках
    plt::plot(freq_vec, ampl_vec, "ro"); // "ro" = червоні кружечки
    plt::show();

	// plt::vlines(freq_vec, std::vector<float>(freq_vec.size(), 0), ampl_vec);
    
    // // додатково можна позначити кружечки на піках
    // plt::scatter(freq_vec, ampl_vec);

    // plt::show();

    // plt::scatter(freq_vec, ampl_vec);
    // plt::show();

#if 0 
	FILE * f_orig = fopen("log_orig.txt", "w"); 
	FILE * f_ft = fopen("log_ft.txt", "w"); 
	FILE * f_mag = fopen("log_mag.txt", "w"); 
	FILE * f_ift = fopen("log_ift.txt", "w"); 
	if(f_orig == NULL)
		perror("Couldn't open file!\n"); 
	if(f_ft == NULL)
		perror("Couldn't open file!\n"); 
	if(f_ift == NULL)
		perror("Couldn't open file!\n"); 
	for(int i = 0; i < N; i++)
		fprintf(f_orig, "%f\n", sound_data_ch0_fl[ptr_shift+i]); 
	for(int i = 0; i < N; i++)
		fprintf(f_mag, "%f\n", magnitudes[i]); 
	// for(int i = 0; i < N; i++)
	// 	fprintf(f_ft, "%f\n", dft_arr[i].real); 
	for(int i = 0; i < N; i++)
		fprintf(f_ift, "%f\n", idft_arr[i].real);  
	fflush(f_orig); 
	fflush(f_ft); 
	fflush(f_mag); 
	fflush(f_ift); 
	fclose(f_orig); 
	fclose(f_ft); 
	fclose(f_mag); 
	fclose(f_ift); 
#endif 

#if 0 
    SDL_Window* window = NULL;
    SDL_Renderer* renderer = NULL;
    SDL_Texture* texture = NULL;
    int i;

    // 1. Initialize SDL
    if (SDL_Init(SDL_INIT_EVERYTHING) < 0) {
        fprintf(stderr, "SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
        return 1;
    }

    // 2. Create window and renderer
    SDL_CreateWindowAndRenderer(SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN, &window, &renderer);
    if (window == NULL || renderer == NULL) {
        fprintf(stderr, "Window or renderer could not be created! SDL_Error: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    // 3. Create a streaming texture
    texture = SDL_CreateTexture(renderer,
                                SDL_PIXELFORMAT_RGBA8888, // Use a common pixel format like RGBA
                                SDL_TEXTUREACCESS_STREAMING,
                                SCREEN_WIDTH, SCREEN_HEIGHT);
    if (texture == NULL) {
        fprintf(stderr, "Texture could not be created! SDL_Error: %s\n", SDL_GetError());
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    // --- Fill your pixel_array with data ---
    // Example: Fill with a gradient or solid color
    for (i = 0; i < SCREEN_WIDTH * SCREEN_HEIGHT; i++) {
        // Example: Red channel varies with x position, green with y
        Uint8 r = (i % SCREEN_WIDTH) * 255 / SCREEN_WIDTH;
        Uint8 g = (i / SCREEN_WIDTH) * 255 / SCREEN_HEIGHT;
        Uint8 b = 100;
        Uint8 a = 255; // Fully opaque

        // Pack the RGBA values into a single Uint32 (adjust order based on your format)
        // For SDL_PIXELFORMAT_RGBA8888, the order in memory might be different due to endianness.
        // The SDL_MapRGBA function helps handle this.
        pixel_array[i] = (a << 24) | (b << 16) | (g << 8) | r;
    }

	// memcpy(pixel_array, sound_data, SCREEN_WIDTH * SCREEN_HEIGHT * sizeof(Uint32)); 

	while (1)
	{ 
        SDL_Event event;
		SDL_PollEvent(&event); 

		if (event.type == SDL_KEYUP)
		{ 
			if (event.key.keysym.sym == SDLK_ESCAPE) 
				break; 

			if (event.quit.type) 
				break; 
		}
			
	
		// 4. Lock the texture to access its internal pixel buffer
		void* mPixels;
		int pitch;
		SDL_LockTexture(texture, NULL, &mPixels, &pitch);

		// 5. Copy your pixel array data to the texture's pixel buffer
		// Ensure the data size matches the texture pitch/dimensions
		memcpy(mPixels, pixel_array, SCREEN_WIDTH * SCREEN_HEIGHT * sizeof(Uint32));

		// 6. Unlock the texture
		SDL_UnlockTexture(texture);
		mPixels = NULL;

		// 7. Render the texture
		SDL_RenderClear(renderer); // Clear the renderer
		SDL_RenderCopy(renderer, texture, NULL, NULL); // Copy the texture to the renderer
		SDL_RenderPresent(renderer); // Update the screen

		// 8. Keep window open for a few seconds
		SDL_Delay(50);
	}

    // 9. Clean up
    SDL_DestroyTexture(texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
#endif 

	free_nc(sound_data); 

	return 0;  
}