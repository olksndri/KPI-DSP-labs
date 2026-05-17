#include "dct.h"
#include <cmath>
#include <math.h>

void dct_1d(float *out, float *inp, int N)
{
    for (int k = 0; k < N; k++)
    {
        float Xk = 0.0f;

        float ak = (k == 0) ? sqrtf(1.0f / N) : sqrtf(2.0f / N);

        for (int n = 0; n < N; n++)
        {
            Xk += inp[n] * cosf((M_PI/N)*(n + 0.5f)*k);
        }

        out[k] = Xk * ak;
    }
}

void idct_1d(float *out, float *inp, int N)
{
	for (int n = 0; n < N; n++)
    {
        float Xn = 0.0f;

        for (int k = 0; k < N; k++)
		{
			float ak = (k == 0) ? sqrtf(1.0f / N) : sqrtf(2.0f / N);

            Xn += ak * inp[k] * cosf((M_PI/N)*(n + 0.5f)*k);
        }

        out[n] = Xn;
    }
}

void dct_2d(float *out, float *inp,  int out_h, int out_w)
{
	int M = out_h;
	int N = out_w;

	for(int v = 0; v < M; v++)
	{
		float Av = (v == 0) ? sqrtf(1.0f / M) : sqrtf(2.0f / M);

		for(int u = 0; u < N; u++)
		{
			float Au = (u == 0) ? sqrtf(1.0f / N) : sqrtf(2.0f / N);

			float acc = 0.0f;

			for(int y = 0; y < M; y++)
			{
				for(int x = 0; x < N; x++)
				{
					float fxy_inp = inp[y * N + x];
					acc += fxy_inp  * cosf( (M_PI/N) * (x + 0.5f) * u)
									* cosf( (M_PI/M) * (y + 0.5f) * v);
				}
			}

			out[v * N + u] = Av * Au * acc;
		}
	}
}

void idct_2d(float *out, float *inp,  int out_h, int out_w)
{
	int M = out_h;
	int N = out_w;

	for(int y = 0; y < M; y++)
	{
		for(int x = 0; x < N; x++)
		{
			float acc = 0.0f;

			for(int v = 0; v < M; v++)
			{
				float Av = (v == 0) ? sqrtf(1.0f / M) : sqrtf(2.0f / M);

				for(int u = 0; u < N; u++)
				{
					float Au = (u == 0) ? sqrtf(1.0f / N) : sqrtf(2.0f / N);

					float fuv_inp = inp[v * N + u];
					acc += Av * Au * fuv_inp  * cosf( (M_PI/N) * (x + 0.5f) * u)
											  * cosf( (M_PI/M) * (y + 0.5f) * v);
				}
			}

			out[y * N + x] = acc;
		}
	}
}
