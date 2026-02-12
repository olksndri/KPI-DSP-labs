#pragma once 

#include <math.h> 

typedef enum WINDOW_T  
{ 
	NO_WINDOW,
	HANN_WINDOW, 
	HAMMING_WINDOW, 
} WINDOW_T; 

void hann_window(float *data, float *windowed_data, int N); 
