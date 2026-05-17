#include <iostream>
#include <vector>
#include <stdint.h>
#include "../dsp/dct.h"
#include "../utils/memops.h"
#include <chrono>
#include <algorithm>
#include <cstdint>
#include <math.h>
#include <stdint.h>
#include <string.h>

typedef enum PAD_TYPE
{
	NO_PAD,
	ZERO,		// Pad with zero value
	REPLICATE, 	// Replicate last value
	REFLECT, 	// Reflects edge pixels
}PAD_TYPE;

typedef struct PADDING
{
	PAD_TYPE padding_type;

	uint32_t top;
	uint32_t bottom;
	uint32_t left;
	uint32_t right;
}PADDING;

typedef enum JPEG_CH_TYPE
{
	LUMINANCE,
	CHROMA,
}JPEG_CH_TYPE;

typedef struct RLE_IMAGE_1CH_EOB
{
	std::vector<int16_t>data;
	std::vector <int8_t>EOB; //last non-zero element index
}RLE_IMAGE_1CH_EOB;



void jpeg_encode(RLE_IMAGE_1CH_EOB &Y_rle, RLE_IMAGE_1CH_EOB &Cb_rle, RLE_IMAGE_1CH_EOB &Cr_rle,
	PADDING &pad_luminance, PADDING &pad_chroma, uint8_t *RGB, int h, int w,  int use_chroma_downsampling);

void jpeg_decode(uint8_t *RGB, RLE_IMAGE_1CH_EOB &Y_rle, RLE_IMAGE_1CH_EOB &Cb_rle, RLE_IMAGE_1CH_EOB &Cr_rle,
	PADDING &pad_luminance, PADDING &pad_chroma, int h, int w,  int use_chroma_downsampling);
