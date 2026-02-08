#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>

#if 1 
#include <SDL2/SDL.h> // for initializing and shutdown functions
#include <SDL2/SDL_image.h> // for rendering images and graphics on screen
#include <SDL2/SDL_timer.h> // for using SDL_Delay() functions
#endif 

#include "include/t.h"


// #include "gnuplot_i.h" // Include the library header





// Розробка програмного забезпечення для швидкого
// перетворення Фур’є голосового сигналу.

// Завдання: розробити програмне забезпечення для
// зчитування аудіоданних із звукового 
// файлу, їх фрагментації та швидкого перетворення Фур’є кожного із
// фрагментів.

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


// #include <SDL2/SDL.h>
// #include <stdio.h>

#define SCREEN_WIDTH 640
#define SCREEN_HEIGHT 480

// Your raw pixel data array (example: 1D array for simplicity)
// The size should be width * height * bytes_per_pixel
// For RGBA (4 bytes per pixel), this would be SCREEN_WIDTH * SCREEN_HEIGHT * 4
Uint32 pixel_array[SCREEN_WIDTH * SCREEN_HEIGHT];

int main()
{
	// printf("hello world ! 1234\n"); 
	// f(); 

	char sf_path[] = "data/ff-16b-2c-44100hz.wav"; 
	
	SNDFILE* f_sound = sf_open(sf_path, SFM_READ, &sf_info);

	log_sf_info(&sf_info); 


	float *sound_data = (float*)malloc_nc(sf_info.frames * sizeof(float)); 

	// sf_count_t read_data = 
	sf_read_float(f_sound, sound_data, sf_info.frames) ;
	// printf("read data: %d\n", read_data); 




#if 1 
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

	memcpy(pixel_array, sound_data, SCREEN_WIDTH * SCREEN_HEIGHT * sizeof(Uint32)); 

	

	int close = 0; 
	while (!close)
	{ 
        SDL_Event event;
		SDL_PollEvent(&event); 

		if (event.type == SDL_KEYDOWN)
			if (event.key.keysym.sym == SDLK_ESCAPE) 
					close = 1;
	
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