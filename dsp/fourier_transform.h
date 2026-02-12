#pragma once 

#include <math.h>

#include "../math/complex_math.h"

void dft(float *in, complex *out, int N); 

void idft(complex *in, complex *out, int N); 