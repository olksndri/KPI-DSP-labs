#pragma once 

typedef struct complex {
	float real; 
	float imag; 
} complex; 

void float_to_complex(float *df, complex *dc, int length); 

complex complex_sub(complex a, complex b);

complex complex_add(complex a, complex b);

complex complex_mul(complex a, complex b); 