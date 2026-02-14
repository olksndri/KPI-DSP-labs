#include "window_funcs.h"


void hann_window_complex(complex *data, complex *windowed_data, int N)
{ 
	assert(N > 1); 

	for(int n = 0; n < N; n++)
	{
		float window =  0.5 * (1 - cosf((2*M_PI*n)/(N-1))); 
		windowed_data[n].real = data[n].real * window; 
		windowed_data[n].imag = data[n].imag * window; 
	}
}

void hann_window(float *data, float *windowed_data, int N)
{ 
	assert(N > 1); 

	for(int n = 0; n < N; n++)
	{
		float window = 0.5 * (1 - cosf((2*M_PI*n)/(N-1))); 
		windowed_data[n] = data[n] * window; 
	}
}

