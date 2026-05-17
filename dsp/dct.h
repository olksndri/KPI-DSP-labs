#pragma once

#include <math.h>
#include <assert.h>

#include "../math/complex_math.h"
#include "../utils/memops.h"


void dct_1d(float *out, float *inp, int N);
void idct_1d(float *out, float *inp, int N);
void dct_2d(float *out, float *inp,  int out_h, int out_w);
void idct_2d(float *out, float *inp,  int out_h, int out_w);
