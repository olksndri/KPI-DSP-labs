#include "fourier_transform.h"


void fft_recursive(complex *in, complex *out, int N) 
{
	assert(N > 0); 

	if(N == 1) 
	{ 
		out[0] = in[0];
		return;
	}

	complex *dft_even = (complex*)malloc_nc(N/2*sizeof(complex)); 
	complex *dft_odd = (complex*)malloc_nc(N/2*sizeof(complex)); 

	int t = 0; 
	for(int n = 0; n < N/2; n++)
	{		
		dft_even[n] = in[n*2]; 
		dft_odd[n] =  in[n*2+1]; 
	} 

	// Recursively calculate both parts of FFT
	fft_recursive(dft_even, dft_even, N/2); 
	fft_recursive(dft_odd, dft_odd, N/2); 


	// Butterfly
	for(int k = 0; k < N / 2; k++)
	{ 
		complex twiddle_factor = 
		{ 
			.real = cosf(2*M_PI*k/N), 
			.imag = -sinf(2*M_PI*k/N), 
		}; 

		out[k] = complex_add(dft_even[k], complex_mul(twiddle_factor, dft_odd[k]));   
		out[k+N/2] = complex_sub(dft_even[k], complex_mul(twiddle_factor, dft_odd[k])); 
	}

	free_nc(dft_even); 
	free_nc(dft_odd); 
} 

void ifft_recursive(complex *in, complex *out, int N) 
{
	assert(N > 0); 

	if(N == 1) 
	{ 
		out[0] = in[0];
		return;
	}

	complex *dft_even = (complex*)malloc_nc(N/2*sizeof(complex)); 
	complex *dft_odd = (complex*)malloc_nc(N/2*sizeof(complex)); 

	int t = 0; 
	for(int n = 0; n < N/2; n++)
	{		
		dft_even[n] = in[n*2]; 
		dft_odd[n] =  in[n*2+1]; 
	} 

	ifft_recursive(dft_even, dft_even, N/2); 
	ifft_recursive(dft_odd, dft_odd, N/2); 

	for(int k = 0; k < N / 2; k++)
	{ 
		complex twiddle_factor = 
		{ 
			.real = cosf(2*M_PI*k/N), 
			.imag = sinf(2*M_PI*k/N), 
		}; 

		out[k] = complex_add(dft_even[k], complex_mul(twiddle_factor, dft_odd[k]));   
		out[k+N/2] = complex_sub(dft_even[k], complex_mul(twiddle_factor, dft_odd[k])); 
		
	}

	free_nc(dft_even); 
	free_nc(dft_odd); 
} 

void ifft_recursive_scale(complex *inout, int N) 
{
	assert(N > 0); 

	float scale = 1.0f / N; 

	for(int n = 0; n < N; n++)
	{ 
		inout[n].real *= scale;
		inout[n].imag *= scale;
	}
} 

void dft(complex *in, complex *out, int N) 
{
	assert(N > 0); 
	
	if(N == 1) 
	{ 
		out[0] = in[0];
		return;
	}

	for(int k = 0; k < N; k++)
	{
		complex Xout = { 0, 0 }; 

		for(int n = 0; n < N; n++)
		{
			complex Xin = in[n];

			complex basis_func = { 
				.real = cosf((2*M_PI*n*k)/N),
				.imag = -sinf((2*M_PI*n*k)/N), 
			}; 

			Xout = complex_add(Xout, complex_mul(Xin, basis_func)); 
		}
		
		out[k] = Xout; 
	}
}

void idft(complex *in, complex *out, int N) 
{
	assert(N > 0); 

	if(N == 1) 
	{ 
		out[0] = in[0];
		return;
	}
	
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

		Xout.real *= 1.0f/N; 
		Xout.imag *= 1.0f/N; 

		out[k] = Xout; 
	}
}
