#pragma once 

#include <math.h> 
#include <assert.h> 

#include "../math/complex_math.h"

typedef enum WINDOW_T  
{ 
	NO_WINDOW,
	HANN_WINDOW, 
	HAMMING_WINDOW, 
} WINDOW_T; 

void hann_window_complex(complex *data, complex *windowed_data, int N); 

void hann_window(float *data, float *windowed_data, int N); 
