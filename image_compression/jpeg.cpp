#include "jpeg.h"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <vector>

void RGBA_to_YCbCr(uint8_t *YCbCr, uint8_t *A, uint8_t *RGBA, int h, int w)
{

	uint8_t *p_out = YCbCr;
	uint8_t *p_alpha = A;
	uint8_t *p_inp = RGBA;

	for(int i = 0; i < h * w; i++)
	{
		float R = p_inp[i * 4 + 0];
		float G = p_inp[i * 4 + 1];
		float B = p_inp[i * 4 + 2];
		float A = p_inp[i * 4 + 3];

		uint8_t Y	= 16 + (65.481 * R + 128.553 * G + 24.966 * B) / 256;
		uint8_t Cb	= 128 + (-37.797 * R - 74.203 * G + 112 * B) / 256;
		uint8_t Cr 	= 128 + (112 * R - 93.786 * G - 18.214 * B) / 256;

		p_out[i * 3 + 0] = Y;
		p_out[i * 3 + 1] = Cb;
		p_out[i * 3 + 2] = Cr;

		p_alpha[i] = A;
	}
}


void YCbCr_to_RGBA(uint8_t *RGBA, uint8_t *YCbCr, uint8_t *A, int h, int w)
{

	uint8_t *p_out = RGBA;
	uint8_t *p_alpha = A;
	uint8_t *p_inp = YCbCr;

	for(int i = 0; i < h * w; i++)
	{
		float Y 	= p_inp[i * 3 + 0];
		float Cb 	= p_inp[i * 3 + 1];
		float Cr	= p_inp[i * 3 + 2];
		uint8_t A 	= p_alpha[i];

		uint8_t R = 298.082 * (Y - 16) / 256.0 + 408.583 * (Cr - 128.0) / 256.0;
		uint8_t G = 298.082 * (Y - 16) / 256.0 - 100.291 * (Cb - 128.0) / 256.0 - 208.120 * (Cr - 128) / 256.0;
		uint8_t B = 298.082 * (Y - 16) / 256.0 + 516.411 * (Cb - 128.0) / 256.0;

		p_out[i * 4 + 0] = R;
		p_out[i * 4 + 1] = G;
		p_out[i * 4 + 2] = B;
		p_out[i * 4 + 3] = A;
	}
}


static inline uint8_t clamp_fl_to_uint8(float v)
{
	// if(v < 0.0f)
	// {
	// 	return 0;
	// }
	// else if (v > UINT8_MAX)
	// {
	// 	return UINT8_MAX;
	// }
	// else
	// {
	// 	return v;
	// }

 int i = (int)roundf(v); // або (int)(v + 0.5f)
    if (i < 0) return 0;
    if (i > 255) return 255;
    return (uint8_t)i;
}
// static inline uint8_t clamp_fl_to_uint8(float v) {

// }

void RGB_to_YCbCr(uint8_t *YCbCr, uint8_t *RGB, int h, int w)
{
	for(int i = 0; i < h * w; i++)
	{
		float R = RGB[i * 3 + 0];
		float G = RGB[i * 3 + 1];
		float B = RGB[i * 3 + 2];

		uint8_t Y	= clamp_fl_to_uint8(16 + (65.481 * R + 128.553 * G + 24.966 * B) / 256.0) ;
		uint8_t Cb	= clamp_fl_to_uint8(128 + (-37.797 * R - 74.203 * G + 112 * B) / 256.0);
		uint8_t Cr 	= clamp_fl_to_uint8(128 + (112 * R - 93.786 * G - 18.214 * B) / 256.0);

		YCbCr[i * 3 + 0] = Y;
		YCbCr[i * 3 + 1] = Cb;
		YCbCr[i * 3 + 2] = Cr;
	}
}


void YCbCr_to_RGB(uint8_t *RGB, uint8_t *YCbCr, int h, int w)
{
	for(int i = 0; i < h * w; i++)
	{
		float Y 	= YCbCr[i * 3 + 0];
		float Cb 	= YCbCr[i * 3 + 1];
		float Cr	= YCbCr[i * 3 + 2];

		uint8_t R = clamp_fl_to_uint8(298.082 * (Y - 16) / 256.0 + 408.583 * (Cr - 128.0) / 256.0);
		uint8_t G = clamp_fl_to_uint8(298.082 * (Y - 16) / 256.0 - 100.291 * (Cb - 128.0) / 256.0 - 208.120 * (Cr - 128) / 256.0);
		uint8_t B = clamp_fl_to_uint8(298.082 * (Y - 16) / 256.0 + 516.411 * (Cb - 128.0) / 256.0);

		RGB[i * 3 + 0] = R;
		RGB[i * 3 + 1] = G;
		RGB[i * 3 + 2] = B;
	}
}

void split_YCbCr(uint8_t *Y, uint8_t *Cb, uint8_t *Cr, const uint8_t *YCbCr, const int h, const int w)
{
	for(int i = 0; i < h * w; i++)
	{
		Y[i] = YCbCr[i*3+0];
		Cb[i] = YCbCr[i*3+1];
		Cr[i] = YCbCr[i*3+2];
	}
}


void concat_YCbCr(uint8_t *YCbCr, const uint8_t *Y, const uint8_t *Cb, const uint8_t *Cr, const int h, const int w)
{
	for(int i = 0; i < h * w; i++)
	{
		YCbCr[i*3+0] = Y[i];
		YCbCr[i*3+1] = Cb[i];
		YCbCr[i*3+2] = Cr[i];
	}
}


// out_h = ceilf(inp_h / 2.0f);
// out_w = ceilf(inp_w / 2.0f);
void downsample_chroma_4_2_0(uint8_t *chroma_ds, const uint8_t *chroma, int inp_h, int inp_w, int out_h, int out_w)
{
	int out_idx = 0;

	for(int y = 0; y < inp_h; y += 2)
	{
		int inp_curr_y = y + 0;
		int inp_next_y = y + 1;

		for(int x = 0; x < inp_w; x += 2)
		{
			float real_pixels = 1.0f;

			int inp_curr_x = x + 0;
			int inp_next_x = x + 1;

			int vtr_flag = (inp_next_x == inp_w);
			int vbl_flag = (inp_next_y == inp_h);
			int vbr_flag = (inp_next_y == inp_h) || (inp_next_x == inp_w);

			real_pixels = vtr_flag ? real_pixels : real_pixels + 1.0f;
			real_pixels = vbl_flag ? real_pixels : real_pixels + 1.0f;
			real_pixels = vbr_flag ? real_pixels : real_pixels + 1.0f;

			float v_top_left  = chroma[inp_curr_y * inp_w + inp_curr_x];
			float v_top_right = vtr_flag ? 0 : chroma[inp_curr_y * inp_w + inp_next_x];
			float v_bot_left  = vbl_flag ? 0 : chroma[inp_next_y * inp_w + inp_curr_x];
			float v_bot_right = vbr_flag ? 0 : chroma[inp_next_y * inp_w + inp_next_x];

			uint8_t v_out = clamp_fl_to_uint8((v_top_left + v_top_right + v_bot_left + v_bot_right) / real_pixels);

			chroma_ds[out_idx] = v_out;
			out_idx++;
		}
	}
}

// should be passed original
// int out_h, int out_w
void upsample_chroma_4_2_0(uint8_t *chroma_us, const uint8_t *chroma, int inp_h, int inp_w, int out_h, int out_w)
{
	for(int y = 0; y < inp_h; y++)
	{
		int out_curr_y = y * 2 + 0;
		int out_next_y = y * 2 + 1;

		for(int x = 0; x < inp_w; x++)
		{
			uint8_t v = chroma[y * inp_w + x];

			int out_curr_x = x * 2 + 0;
			int out_next_x = x * 2 + 1;

			int out_tl_idx = out_curr_y * out_w + out_curr_x;
			int out_tr_idx = out_curr_y * out_w + out_next_x;
			int out_bl_idx = out_next_y * out_w + out_curr_x;
			int out_br_idx = out_next_y * out_w + out_next_x;

			int vtr_flag = (out_next_x == out_w);
			int vbl_flag = (out_next_y == out_h);
			int vbr_flag = (out_next_y == out_h) || (out_next_x == out_w);

			chroma_us[out_tl_idx] = v;
			if(!vtr_flag) chroma_us[out_tr_idx] = v;
			if(!vbl_flag) chroma_us[out_bl_idx] = v;
			if(!vbr_flag) chroma_us[out_br_idx] = v;
		}
	}
}

void uint8_to_int8(int8_t *out, uint8_t *inp, int len)
{
	for(int i = 0; i < len; i++)
	{
		int16_t v = inp[i];
		out[i] = v - 128;
	}
}

void int8_to_uint8(uint8_t *out, int8_t *inp, int len)
{
	for(int i = 0; i < len; i++)
	{
		int16_t v = inp[i];
		out[i] = v + 128;
	}
}

void int8_to_float(float *out, int8_t *inp, int len)
{
	for(int i = 0; i < len; i++)
		out[i] = inp[i];
}


void float_to_int8(int8_t *out, float *inp, int len)
{
	for(int i = 0; i < len; i++)
		out[i] = roundf(inp[i]);
}


// // Matrix for Y (Luminance)
const uint8_t std_lum_quant[64] = {
    16, 11, 10, 16, 24, 40, 51, 61,
    12, 12, 14, 19, 26, 58, 60, 55,
    14, 13, 16, 24, 40, 57, 69, 56,
    14, 17, 22, 29, 51, 87, 80, 62,
    18, 22, 37, 56, 68, 109, 103, 77,
    24, 35, 55, 64, 81, 104, 113, 92,
    49, 64, 78, 87, 103, 121, 120, 101,
    72, 92, 95, 98, 112, 100, 103, 99
};

// // Matrix for Cb, Cr (Chrominance)
const uint8_t std_chrom_quant[64] = {
    17, 18, 24, 47, 99, 99, 99, 99,
    18, 21, 26, 66, 99, 99, 99, 99,
    24, 26, 56, 99, 99, 99, 99, 99,
    47, 66, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99
};

const uint8_t zig_zag_map[64] = {
     0,  1,  8, 16,  9,  2,  3, 10,
    17, 24, 32, 25, 18, 11,  4,  5,
    12, 19, 26, 33, 40, 48, 41, 34,
    27, 20, 13,  6,  7, 14, 21, 28,
    35, 42, 49, 56, 57, 50, 43, 36,
    29, 22, 15, 23, 30, 37, 44, 51,
    58, 59, 52, 45, 38, 31, 39, 46,
    53, 60, 61, 54, 47, 55, 62, 63
};

void apply_zigzag(int16_t *out_vector, int16_t *in_8x8_block) {
	for (int i = 0; i < 64; i++) {
		out_vector[i] = in_8x8_block[zig_zag_map[i]];
	}
}

void deapply_zigzag(int16_t *out_8x8_block, int16_t *inp_vector)
{
	for (int i = 0; i < 64; i++) {
		out_8x8_block[zig_zag_map[i]] = inp_vector[i];
	}
}

void quantize(int16_t *quantized, float *inp, const uint8_t *q_table, int len)
{
	for(int i = 0; i < len; i++)
		quantized[i] = roundf(inp[i] / q_table[i]);
}

void dequantize(float *dequantized, int16_t *inp, const uint8_t *q_table, int len)
{
	for(int i = 0; i < len; i++)
		dequantized[i] = roundf((float)inp[i] * q_table[i]);
}

// Assumes 1 channel image as input
void unpadding(uint8_t *unpadded_img, uint8_t *padded_img, int orig_h, int orig_w, PADDING &pad)
{
	if(pad.padding_type == PAD_TYPE::NO_PAD)
		return;

	int pad_h = pad.top + pad.bottom;
	int pad_w = pad.left + pad.right;

	int padded_h = orig_h + pad_h;
	int padded_w = orig_w + pad_w;

	// Copy input image into padded image array
	for(int y = 0; y < orig_h; y++)
	{
		uint8_t* p_padded_img_row = &padded_img[(pad.top * padded_w + pad.left) + y * padded_w];
		uint8_t* p_upadded_img_row = &unpadded_img[y * orig_w];
		memcpy(p_upadded_img_row, p_padded_img_row, orig_w * sizeof(uint8_t));
	}
}

// Assumes 1 channel image as input
void padding(uint8_t *padded_img, uint8_t *inp_img, int orig_h, int orig_w, PADDING &pad)
{
	if(pad.padding_type == PAD_TYPE::NO_PAD)
		return;

	int pad_h = pad.top + pad.bottom;
	int pad_w = pad.left + pad.right;

	int padded_h = orig_h + pad_h;
	int padded_w = orig_w + pad_w;

	// Copy input image into padded image array
	for(int y = 0; y < orig_h; y++)
	{
		uint8_t* p_padded_img_row = &padded_img[(pad.top * padded_w + pad.left) + y * padded_w];
		uint8_t* p_inp_img_row = &inp_img[y * orig_w];
		memcpy(p_padded_img_row, p_inp_img_row, orig_w * sizeof(uint8_t));
	}


	if(pad.padding_type == PAD_TYPE::REFLECT)
	{
		// Left and Right padding
		for(int y = pad.top; y < (pad.top + orig_h); y++)
		{
			if(pad.left > 0)
			{
				for(int x = 0; x < pad.left; x++)
				{
					int x_inp_idx = (pad.left * 2 - 1) - x;
					int out_idx = y * padded_w + x;
					int inp_idx = y * padded_w + x_inp_idx;
					padded_img[out_idx] = padded_img[inp_idx];
				}
			}

			if(pad.right > 0)
			{
				int x_start = pad.left + orig_w;
				for(int x = x_start; x < padded_w; x++)
				{
					int x_inp_idx = (x_start - 1) - (x - x_start);
					int out_idx = y * padded_w + x;
					int inp_idx = y * padded_w + x_inp_idx;
					padded_img[out_idx] = padded_img[inp_idx];
				}
			}
		}

		// Top padding
		if(pad.top)
		{
			for(int y = 0; y < pad.top; y++)
			{
				int y_inp_idx = (pad.top * 2 - 1) - y;
				uint8_t *p_inp = &padded_img[y_inp_idx * padded_w];
				uint8_t *p_out = &padded_img[y * padded_w];
				memcpy(p_out, p_inp, padded_w * sizeof(uint8_t));
			}
		}

		// Bottom padding
		if(pad.bottom)
		{
			int y_start = pad.top + orig_h;
			for(int y = y_start; y < padded_h; y++)
			{
				int y_inp_idx = (y_start - 1) - (y - y_start);
				uint8_t *p_inp = &padded_img[y_inp_idx * padded_w];
				uint8_t *p_out = &padded_img[y * padded_w];
				memcpy(p_out, p_inp, padded_w * sizeof(uint8_t));
			}
		}
	}
	else
	{
		printf("Not implemented yet!\n");
		exit(0);
	}
}

void jpeg_process_channel_forward(int16_t *out, float *color, int block_h, int block_w, int h, int w, JPEG_CH_TYPE ch_type)
{
	std::vector<float> block(block_h * block_w);
	float *p_block = block.data();

	std::vector<float> block_dct2d(block_h * block_w);
	std::vector<int16_t> block_quantized(block_h * block_w);

	const uint8_t *q_table = (ch_type == JPEG_CH_TYPE::LUMINANCE) ? std_lum_quant :  std_chrom_quant;

	int16_t *p_temp_out = out;
	for(int y = 0; y < h; y += block_h)
	{
		for(int x = 0; x < w; x += block_w)
		{
			// Copy block_h x block_w block from input image
			for(int by = 0; by < block_w; by++)
				for(int bx = 0; bx < block_h; bx++)
					p_block[by * block_w + bx] = color[(y + by) * w + x + bx];

			// Perform forward 2D Discrete Cosine Transform to get frequency coefficients of signal
			// In result: first value is a DC component - avarage brightness of image
			// and other array values are freq. components from lower to higher frequency
			// The further to the right and lower, the higher the frequency (fine details, sharp edges)
			// that represented by a coefficient.
		 	dct_2d(block_dct2d.data(), block.data(), block_h, block_w);

			// Apply quantization to the frequency coefficient whith such a table
			// that will result into a zeroing all high-freq components, generating
			// a large number of "0"
			quantize(block_quantized.data(), block_dct2d.data(), q_table, block_h * block_w);

			// Apply zigzag scanning to reorder data such way
			// that end of array will be filled with sequence of zeros
			// apply_zigzag(&out[y * w + x], block_quantized.data());
			apply_zigzag(p_temp_out, block_quantized.data());
			p_temp_out+=64;

		}
	}

}

void jpeg_process_channel_backward(float *color, int16_t *compressed_in, int block_h, int block_w, int h, int w, JPEG_CH_TYPE ch_type)
{
	std::vector<float> block(block_h * block_w);
	float *p_block = block.data();

	std::vector<float> block_dct2d(block_h * block_w);
	std::vector<int16_t> block_quantized(block_h * block_w);

	const uint8_t *q_table = (ch_type == JPEG_CH_TYPE::LUMINANCE) ? std_lum_quant :  std_chrom_quant;
	int16_t *p_temp_in = compressed_in;

	for(int y = 0; y < h; y += block_h)
	{
		for(int x = 0; x < w; x += block_w)
		{
			// Apply zigzag scanning to reorder data such way
			// that end of array will be filled with sequence of zeros
			// deapply_zigzag(block_quantized.data(), compressed_in);
			deapply_zigzag(block_quantized.data(), p_temp_in);
			p_temp_in+=64 ;
			// Apply quantization to the frequency coefficient whith such a table
			// that will result into a zeroing all high-freq components, generating
			// a large number of "0"
			dequantize(block_dct2d.data(), block_quantized.data(), q_table, block_h * block_w);

			// Perform forward 2D Discrete Cosine Transform to get frequency coefficients of signal
			// In result: first value is a DC component - avarage brightness of image
			// and other array values are freq. components from lower to higher frequency
			// The further to the right and lower, the higher the frequency (fine details, sharp edges)
			// that represented by a coefficient.
		 	idct_2d(block.data(), block_dct2d.data(), block_h, block_w);

			// Copy block_h x block_w block from input image
			for(int by = 0; by < block_w; by++)
				for(int bx = 0; bx < block_h; bx++)
					color[(y + by) * w + x + bx] = p_block[by * block_w + bx];
		}
	}
}



void rle_compress_eob(RLE_IMAGE_1CH_EOB &out, int16_t *inp, int block_h, int block_w, int h, int w)
{
	int block_size = block_h * block_w;
	int blocks_count = (h * w) / block_size;

	for(int i = 0; i < blocks_count; i++)
	{
		int start_idx = ((i + 1) * block_size) - 1;
		int end_idx = i * block_size;

		int eob_pushed = 0;
		int eob_cnt = block_size - 1;
		for(int b = start_idx; b >= end_idx; b--)
		{
			if(inp[b] != 0)
			{
				int eob = eob_cnt;
				out.EOB.push_back(eob);
				for(int nz = 0; nz <= eob; nz++)
				{
					out.data.push_back(inp[i * block_size + nz]);
				}
				eob_pushed = 1;
				break;
			}

			eob_cnt--;
		}

		// special value to specify that all the block contains zeros at end
		if(!eob_pushed)
		{
			out.EOB.push_back(-1);
		}
	}
}


void rle_decompress_eob(int16_t *out, RLE_IMAGE_1CH_EOB &inp, int block_h, int block_w, int h, int w)
{
	int block_size = block_h * block_w;
	int blocks_count = (h * w) / block_size;
	int cur_inp_block_idx = 0;
	for(int i = 0; i < blocks_count; i++)
	{
		int cur_eob = inp.EOB[i];

		int cur_out_block_idx = i * block_size;
		if(cur_eob == -1)
		{
			memset(&out[cur_out_block_idx], 0x0, block_size * sizeof(int16_t));
		}
		else if(cur_eob == block_size - 1)
		{
			memcpy(&out[cur_out_block_idx], &inp.data[cur_inp_block_idx], block_size * sizeof(int16_t));
		}
		else // block contains zero (s) at the end and non-zeros at the begin
		{
			int non_zero_elements = (cur_eob + 1);
			int zero_elements = block_size - non_zero_elements;
			memcpy(&out[cur_out_block_idx], &inp.data[cur_inp_block_idx], non_zero_elements * sizeof(int16_t));
			memset(&out[cur_out_block_idx + non_zero_elements], 0x0, zero_elements * sizeof(int16_t));
		}

		// if cur_eob == -1 will be incremented by zero because 0 elements are non-zero
		// for current block
		cur_inp_block_idx += cur_eob+1;
	}
}

void print_img_numbers(const char *fname, int16_t *t, int h, int w)
{
	FILE *fp = fopen(fname, "w");
	for(int y = 0; y < h; y++)
	{
		fprintf(fp, "\n");
		for(int x = 0; x < w; x++)
		{
			fprintf(fp, "%d ", t[y*w+x]);
		}
	}
	fflush(fp);
	fclose(fp);
}

void jpeg_encode(RLE_IMAGE_1CH_EOB &Y_rle, RLE_IMAGE_1CH_EOB &Cb_rle, RLE_IMAGE_1CH_EOB &Cr_rle,
	PADDING &pad_luminance, PADDING &pad_chroma, uint8_t *RGB, int h, int w, int use_chroma_downsampling)
{
	std::vector<uint8_t> YCbCr(h * w * 3);

	// The image is converted from RGB to YCbCr
	RGB_to_YCbCr(YCbCr.data(), RGB, h, w);

	// Split YCbCr image into three separate images
	std::vector<uint8_t> Y(h * w);
	std::vector<uint8_t> Cb(h * w);
	std::vector<uint8_t> Cr(h * w);
	split_YCbCr(Y.data(), Cb.data(), Cr.data(), YCbCr.data(), h, w);

	// Chroma channels could be downsampled, resulting into 2 times lesser height and width
	int h_ds = (h + 1) / 2;
	int w_ds = (w + 1) / 2;
	std::vector<uint8_t> Cb_ds(h_ds * w_ds);
	std::vector<uint8_t> Cr_ds(h_ds * w_ds);

	if(use_chroma_downsampling)
	{
		downsample_chroma_4_2_0(Cb_ds.data(), Cb.data(), h, w, h_ds, w_ds);
		downsample_chroma_4_2_0(Cr_ds.data(), Cr.data(), h, w, h_ds, w_ds);
	}
	else
	{
		h_ds = h;
		w_ds = w;
		Cb_ds.resize(h_ds * w_ds);
		Cr_ds.resize(h_ds * w_ds);
		std::copy(Cb.begin(), Cb.end(), Cb_ds.begin());
		std::copy(Cr.begin(), Cr.end(), Cr_ds.begin());
	}
	// Pad images to have a multiple of 8 size
	int mulof8_h = (h + 7) & -8;
	int mulof8_w = (w + 7) & -8;
	int mulof8_hds = (h_ds + 7) & -8;
	int mulof8_wds = (w_ds + 7) & -8;

	std::vector<uint8_t> Y_padded(mulof8_h * mulof8_w);
	std::vector<uint8_t> Cb_padded(mulof8_hds * mulof8_wds);
	std::vector<uint8_t> Cr_padded(mulof8_hds * mulof8_wds);

	if((mulof8_h != h) || (mulof8_w != w))
	{
		uint32_t padr = mulof8_w - w;
		uint32_t padb = mulof8_h - h;
		PADDING pad = { .padding_type = PAD_TYPE::REFLECT, .top = 0, .bottom = padb, .left = 0, .right = padr };
		pad_luminance = pad;
		padding(Y_padded.data(), Y.data(), h, w, pad);
	}
	else
	{
		std::copy(Y.begin(), Y.end(), Y_padded.begin());
		pad_luminance = { .padding_type = PAD_TYPE::NO_PAD, .top = 0, .bottom = 0, .left = 0, .right = 0} ;
	}
	if((mulof8_hds != h_ds) || (mulof8_wds != w_ds))
	{
		uint32_t padr = mulof8_wds - w_ds;
		uint32_t padb = mulof8_hds - h_ds;
		PADDING pad = { .padding_type = PAD_TYPE::REFLECT, .top = 0, .bottom = padb, .left = 0, .right = padr };
		pad_chroma = pad;
		padding(Cb_padded.data(), Cb_ds.data(), h_ds, w_ds, pad);
		padding(Cr_padded.data(), Cr_ds.data(), h_ds, w_ds, pad);
	}
	else
	{
		std::copy(Cb_ds.begin(), Cb_ds.end(), Cb_padded.begin());
		std::copy(Cr_ds.begin(), Cr_ds.end(), Cr_padded.begin());
		pad_chroma = { .padding_type = PAD_TYPE::NO_PAD, .top = 0, .bottom = 0, .left = 0, .right = 0} ;
	}

	std::vector<int8_t> Y_pad_int8(Y_padded.size());
	std::vector<int8_t> Cb_pad_int8(Cb_padded.size());
	std::vector<int8_t> Cr_pad_int8(Cr_padded.size());
	uint8_to_int8(Y_pad_int8.data(), Y_padded.data(), Y_padded.size());
	uint8_to_int8(Cb_pad_int8.data(), Cb_padded.data(), Cb_padded.size());
	uint8_to_int8(Cr_pad_int8.data(), Cr_padded.data(), Cr_padded.size());

	std::vector<float> Y_pad_fl(Y_padded.size());
	std::vector<float> Cb_pad_fl(Cb_padded.size());
	std::vector<float> Cr_pad_fl(Cr_padded.size());
	int8_to_float(Y_pad_fl.data(), Y_pad_int8.data(), Y_pad_int8.size());
	int8_to_float(Cb_pad_fl.data(), Cb_pad_int8.data(), Cb_pad_int8.size());
	int8_to_float(Cr_pad_fl.data(), Cr_pad_int8.data(), Cr_pad_int8.size());

	// Break image channels into 8x8 blocks and process each block separately.
	int block_h = 8;
	int block_w = 8;
 	std::vector<int16_t> Y_zg(Y_padded.size());
 	std::vector<int16_t> Cb_zg(Cb_padded.size());
 	std::vector<int16_t> Cr_zg(Cr_padded.size());
	jpeg_process_channel_forward(Y_zg.data(), Y_pad_fl.data(), block_h, block_w, mulof8_h, mulof8_w, JPEG_CH_TYPE::LUMINANCE);
	jpeg_process_channel_forward(Cb_zg.data(), Cb_pad_fl.data(), block_h, block_w, mulof8_hds, mulof8_wds, JPEG_CH_TYPE::CHROMA);
	jpeg_process_channel_forward(Cr_zg.data(), Cr_pad_fl.data(), block_h, block_w, mulof8_hds, mulof8_wds, JPEG_CH_TYPE::CHROMA);

	rle_compress_eob(Y_rle, Y_zg.data(), block_h, block_w, mulof8_h, mulof8_w);
	rle_compress_eob(Cb_rle, Cb_zg.data(), block_h, block_w, mulof8_hds, mulof8_wds);
	rle_compress_eob(Cr_rle, Cr_zg.data(), block_h, block_w, mulof8_hds, mulof8_wds);
}

void jpeg_decode(uint8_t *RGB, RLE_IMAGE_1CH_EOB &Y_rle, RLE_IMAGE_1CH_EOB &Cb_rle, RLE_IMAGE_1CH_EOB &Cr_rle,
	PADDING &pad_luminance, PADDING &pad_chroma, int h, int w, int use_chroma_downsampling)
{
	// Start restoring image by concatenating blocks back to image channels
	int h_ds = (h + 1) / 2;
	int w_ds = (w + 1) / 2;

	if(!use_chroma_downsampling)
	{
		h_ds = h;
		w_ds = w;
	}

	// Pad images to have a multiple of 8 size
	int mulof8_h = (h + 7) & -8;
	int mulof8_w = (w + 7) & -8;
	int mulof8_hds = (h_ds + 7) & -8;
	int mulof8_wds = (w_ds + 7) & -8;

	std::vector<uint8_t> Y_padded(mulof8_h * mulof8_w);
	std::vector<uint8_t> Cb_padded(mulof8_hds * mulof8_wds);
	std::vector<uint8_t> Cr_padded(mulof8_hds * mulof8_wds);

	int block_h = 8;
	int block_w = 8;

	std::vector<int16_t> Y_zg(mulof8_h * mulof8_w);
	std::vector<int16_t> Cb_zg(mulof8_hds * mulof8_wds);
	std::vector<int16_t> Cr_zg(mulof8_hds * mulof8_wds);

	rle_decompress_eob(Y_zg.data(), Y_rle, block_h, block_w, mulof8_h, mulof8_w);
	rle_decompress_eob(Cb_zg.data(), Cb_rle, block_h, block_w, mulof8_hds, mulof8_wds);
	rle_decompress_eob(Cr_zg.data(), Cr_rle, block_h, block_w, mulof8_hds, mulof8_wds);

	std::vector<float> Y_pad_fl(Y_zg.size());
	std::vector<float> Cb_pad_fl(Cb_zg.size());
	std::vector<float> Cr_pad_fl(Cr_zg.size());
	jpeg_process_channel_backward(Y_pad_fl.data(), Y_zg.data(), block_h, block_w, mulof8_h, mulof8_w, JPEG_CH_TYPE::LUMINANCE);
	jpeg_process_channel_backward(Cb_pad_fl.data(), Cb_zg.data(), block_h, block_w, mulof8_hds, mulof8_wds, JPEG_CH_TYPE::CHROMA);
	jpeg_process_channel_backward(Cr_pad_fl.data(), Cr_zg.data(), block_h, block_w, mulof8_hds, mulof8_wds, JPEG_CH_TYPE::CHROMA);

	std::vector<int8_t> Y_pad_int8(Y_padded.size());
	std::vector<int8_t> Cb_pad_int8(Cb_padded.size());
	std::vector<int8_t> Cr_pad_int8(Cr_padded.size());
	float_to_int8(Y_pad_int8.data(), Y_pad_fl.data(), Y_pad_int8.size());
	float_to_int8(Cb_pad_int8.data(), Cb_pad_fl.data(), Cb_pad_int8.size());
	float_to_int8(Cr_pad_int8.data(), Cr_pad_fl.data(), Cr_pad_int8.size());

	int8_to_uint8(Y_padded.data(), Y_pad_int8.data(), Y_padded.size());
	int8_to_uint8(Cb_padded.data(), Cb_pad_int8.data(), Cb_padded.size());
	int8_to_uint8(Cr_padded.data(), Cr_pad_int8.data(), Cr_padded.size());

	// Split YCbCr image into three separate images
	std::vector<uint8_t> Y(h * w);
	std::vector<uint8_t> Cb(h * w);
	std::vector<uint8_t> Cr(h * w);

	// Chroma channels could be downsampled, resulting into 2 times lesser height and width
	std::vector<uint8_t> Cb_ds(h_ds * w_ds);
	std::vector<uint8_t> Cr_ds(h_ds * w_ds);

	unpadding(Y.data(), Y_padded.data(), h, w, pad_luminance);
	unpadding(Cb_ds.data(), Cb_padded.data(), h_ds, w_ds, pad_chroma);
	unpadding(Cr_ds.data(), Cr_padded.data(), h_ds, w_ds, pad_chroma);

	if(use_chroma_downsampling)
	{
		upsample_chroma_4_2_0(Cb.data(), Cb_ds.data(), h_ds, w_ds, h, w);
		upsample_chroma_4_2_0(Cr.data(), Cr_ds.data(), h_ds, w_ds, h, w);
	}
	else
	{
		std::copy(Cb_ds.begin(), Cb_ds.end(), Cb.begin());
		std::copy(Cr_ds.begin(), Cr_ds.end(), Cr.begin());
	}

	std::vector<uint8_t> YCbCr(h * w * 3);

	concat_YCbCr(YCbCr.data(), Y.data(), Cb.data(), Cr.data(), h, w);

	YCbCr_to_RGB(RGB, YCbCr.data(), h, w);
}
