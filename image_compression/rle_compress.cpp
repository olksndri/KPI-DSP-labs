#include "rle.h"
#include <cstdint>
#include <cstring>

// Assumes that input image has 4 channels (e.g. RGBA or any other 4ch format in any order)
void rle_compress(uint32_t *inp, RLE_IMAGE &st, int w, int h)
{
	RLE_PX px;

	uint32_t times_size_bytes = sizeof(px.times);
	uint32_t times_max;

	if(times_size_bytes == 1)
		times_max = UINT8_MAX;
	else if(times_size_bytes == 2)
		times_max = UINT16_MAX;
	else if(times_size_bytes == 4)
		times_max = UINT32_MAX;
	else
	{
		perror("rle_encode error: no support for RLE_PX.times declared as uint64_t type\n");
		exit(0);
	}

	uint32_t max_out_size_bytes = w * h * sizeof(RLE_PX);
	RLE_PX *temp = (RLE_PX *)malloc(max_out_size_bytes);

	const uint32_t *p_inp = inp;
	const uint32_t *p_inp_max = inp + w * h;
	RLE_PX *p_temp = temp;

	int i = 0;
	for(; i < w * h; )
	{
		uint32_t cur_val = p_inp[0];
		uint32_t seq_times = 1;

		const uint32_t *pp_inp = p_inp+1;
		while(	pp_inp < p_inp_max &&
				seq_times < times_max &&
		 		cur_val == pp_inp[0] )
		{
			pp_inp++;
			seq_times++;
		}

		p_inp = pp_inp;

		px = { .value = cur_val, .times = seq_times };
		*p_temp = px;
		p_temp++;

		i += seq_times;
	}

	int actual_out_size = p_temp - temp;
	uint32_t actual_out_size_bytes = actual_out_size * sizeof(RLE_PX);

	st.data = (RLE_PX *)malloc(actual_out_size_bytes);
	memcpy(st.data, temp, actual_out_size_bytes);

	st.data_size = actual_out_size;
	st.orig_w = w;
	st.orig_h = h;

	free(temp);
}
