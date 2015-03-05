#include "Chunk.h"

#ifndef __COMMON_CUDA_H__
#define __COMMON_CUDA_H__
typedef struct _DataContainer
{
    float* u;
    float* v;
    float* w;
    float* freq;
    float* data_real;
    float* data_imag;
    float* data_weight;
    int* spw;
	int* field;
} DataContainer;

void copy_data_to_cuda(DataContainer& data, Chunk& chunk);
void copy_data_to_host(DataContainer& data, Chunk& chunk);
void setup_freq(DataContainer& data, float* freq,
                const size_t nchan, const size_t nspw);

#endif // inclusion guard
