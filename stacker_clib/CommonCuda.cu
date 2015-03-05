#include "cuda_error.h"
#include "Chunk.h"
#include "CommonCuda.h"
void setup_freq(DataContainer& data, float* freq, const size_t nchan,/*{{{*/
                const size_t nspw)
{
    CudaSafeCall(cudaMalloc( (void**)&data.freq, sizeof(float)*nchan*nspw));
    CudaSafeCall(cudaMemcpy(data.freq, freq, sizeof(float)*nchan*nspw,
                cudaMemcpyHostToDevice));
}/*}}}*/
void copy_data_to_cuda(DataContainer& data, Chunk& chunk)/*{{{*/
{
    size_t chunk_size = chunk.size();

    float* u = new float[chunk_size];
    float* v = new float[chunk_size];
    float* w = new float[chunk_size];
    int* spw = new int[chunk_size];
    int* field = new int[chunk_size];
    for(size_t uvrow = 0; uvrow < chunk_size; uvrow++)
    {
        u[uvrow] = chunk.inVis[uvrow].u;
        v[uvrow] = chunk.inVis[uvrow].v;
        w[uvrow] = chunk.inVis[uvrow].w;
        spw[uvrow] = chunk.inVis[uvrow].spw;
        field[uvrow] = chunk.inVis[uvrow].fieldID;
    }

    CudaSafeCall(cudaMemcpy(data.u, u, sizeof(float)*chunk_size,
                cudaMemcpyHostToDevice));
    CudaSafeCall(cudaMemcpy(data.v, v, sizeof(float)*chunk_size,
                cudaMemcpyHostToDevice));
    CudaSafeCall(cudaMemcpy(data.w, w, sizeof(float)*chunk_size,
                cudaMemcpyHostToDevice));
    CudaSafeCall(cudaMemcpy(data.spw, spw, sizeof(float)*chunk_size,
                cudaMemcpyHostToDevice));
    CudaSafeCall(cudaMemcpy(data.field, field, sizeof(float)*chunk_size,
                cudaMemcpyHostToDevice));
    CudaSafeCall(cudaMemcpy(data.data_real,   chunk.data_real_in,
                sizeof(float)*chunk.size()*chunk.nChan()*chunk.nStokes(),
                cudaMemcpyHostToDevice));
    CudaSafeCall(cudaMemcpy(data.data_imag,   chunk.data_imag_in,
                sizeof(float)*chunk.size()*chunk.nChan()*chunk.nStokes(),
                cudaMemcpyHostToDevice));
    CudaSafeCall(cudaMemcpy(data.data_weight, chunk.weight_in,
                sizeof(float)*chunk.size()*chunk.nStokes(),
                cudaMemcpyHostToDevice));

    delete[] u;
    delete[] v;
    delete[] w;
    delete[] spw;
    delete[] field;
}/*}}}*/
void copy_data_to_host(DataContainer& data, Chunk& chunk)/*{{{*/
{
    CudaSafeCall(cudaMemcpy(chunk.data_real_out, data.data_real,
                sizeof(float)*chunk.size()*chunk.nChan()*chunk.nStokes(),
                cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaMemcpy(chunk.data_imag_out, data.data_imag,
                sizeof(float)*chunk.size()*chunk.nChan()*chunk.nStokes(),
                cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaMemcpy(chunk.weight_out, data.data_weight,
                sizeof(float)*chunk.size()*chunk.nStokes(),
                cudaMemcpyDeviceToHost));
}/*}}}*/
