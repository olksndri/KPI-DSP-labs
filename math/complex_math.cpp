#include "complex_math.h"

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
