#include "fourier_transform.h"

void dft(float *in, complex *out, int N) 
{
	for(int k = 0; k < N; k++)
	{
		complex Xout = { 0, 0 }; 

		for(int n = 0; n < N; n++)
		{
			float Xin = in[n];

			Xout.real += cosf((2*M_PI*n*k)/N) * Xin; 
			Xout.imag += -sinf((2*M_PI*n*k)/N) * Xin; 
		}
		
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
