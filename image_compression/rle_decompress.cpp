#include "rle.h"

// Assumes that input image has 4 channels (e.g. RGBA or any other 4ch format in any order)
void rle_decompress(uint32_t *out, RLE_IMAGE &st)
{
	const RLE_PX *p_inp = st.data;
	uint32_t *p_out = out;

	for(int i = 0; i < st.data_size; i++)
	{
		RLE_PX v = st.data[i];

		for(int k = 0; k < v.times; k++)
		{
			*p_out = v.value;
			p_out++;
		}
	}
}
