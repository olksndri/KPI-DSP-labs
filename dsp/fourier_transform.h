#pragma once 

#include <math.h>
#include <assert.h> 

#include "../math/complex_math.h"
#include "../utils/memops.h"

void fft_recursive(complex *in, complex *out, int N); 

void ifft_recursive(complex *in, complex *out, int N);

void ifft_recursive_scale(complex *inout, int N); 

void dft(complex *in, complex *out, int N); 

void idft(complex *in, complex *out, int N); 