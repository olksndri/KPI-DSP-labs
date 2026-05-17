#include "stdint.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include <vector>


typedef struct RLE_PX
{
	uint32_t value;
	uint8_t times;
}RLE_PX;

typedef struct RLE_IMAGE
{
	RLE_PX * data;
	uint32_t data_size;
	uint32_t orig_w;
	uint32_t orig_h;
}RLE_IMAGE;


void rle_compress(uint32_t *inp, RLE_IMAGE &st, int w, int h);

void rle_decompress(uint32_t *out, RLE_IMAGE &st);
