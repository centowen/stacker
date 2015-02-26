// stacker, Python module for stacking of interferometric data.
// Copyright (C) 2014  Lukas Lindroos
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. 
//
// Library to stack and modsub ms data.
#include "StackChunkComputerGpu_cuda.h"
#include "DataIO.h"
#include "PrimaryBeam.h"
#include "Chunk.h"
#include "Coords.h"
#include "definitions.h"
#include "cuda_error.h"

__constant__ float dev_omega[3*N_MAX_COORDS];
__constant__ float dev_weight[N_MAX_COORDS];

// __global__ void cuVisStack(DataContainer data, CoordContainer coords, 
//                          int chunk_size, int nchan, int n_coords);

__global__ void cuVisStack(DataContainer data, CoordContainer coords, /*{{{*/
                           size_t chunk_size, size_t nchan, size_t nstokes)
{
    size_t uvrow = threadIdx.x + blockIdx.x*blockDim.x;

    float exponent, sin_exponent, cos_exponent, pbcor;
    float m_num_imag, m_num_real, norm;
	float m_real, m_imag;
    float data_real_buff, data_imag_buff;

    while(uvrow < chunk_size)
    {
        float* freq = &data.freq[data.spw[uvrow]*nchan];
        for(size_t chanID = 0; chanID < nchan; chanID++)
        {
            m_num_real = 0.;
            m_num_imag = 0.;
            norm = 0.;

            for(size_t posID = 0; posID < coords.n_coords; posID++)
            {
                size_t pbindex = data.spw[uvrow]*nchan*coords.n_coords
                               + chanID*coords.n_coords
                               + posID;
                pbcor = coords.pb[pbindex];
                exponent = -freq[chanID]*(dev_omega[posID*3+0]*data.u[uvrow] + 
                                          dev_omega[posID*3+1]*data.v[uvrow] + 
                                          dev_omega[posID*3+2]*data.w[uvrow]);
                sincos(exponent, &sin_exponent, &cos_exponent);
                m_num_real += dev_weight[posID]*pbcor*cos_exponent;
                m_num_imag += dev_weight[posID]*pbcor*sin_exponent;
                norm += pbcor*pbcor*dev_weight[posID];
            }

            m_real = m_num_real/norm;
            m_imag = m_num_imag/norm;

            for(size_t stokesID = 0; stokesID < nstokes; stokesID++)
            {
                size_t weightindex = uvrow*nstokes+stokesID;
                size_t dataindex = nchan*weightindex + chanID;

                data_real_buff = data.data_real[dataindex]*m_real -
                                 data.data_imag[dataindex]*m_imag;
                data_imag_buff = data.data_real[dataindex]*m_imag +
                                 data.data_imag[dataindex]*m_real;
//                 data_real_buff = 0.;
//                 data_imag_buff = 0.;
//                 data_real_buff = dev_omega[0];
//                 data_imag_buff = dev_omega[1];

				data.data_real[dataindex] = data_real_buff;
				data.data_imag[dataindex] = data_imag_buff;
            }
        }

        uvrow += blockDim.x*gridDim.x;
    }
};/*}}}*/

void allocate_cuda_data(DataContainer& data, CoordContainer& dev_coords,/*{{{*/
                        const size_t nchan, const size_t nstokes,
                        const size_t chunk_size, const size_t nmaxcoords,
                        const size_t nspw)
{
    CudaSafeCall(cudaMalloc( (void**)&data.u, sizeof(float)*chunk_size));
    CudaSafeCall(cudaMalloc( (void**)&data.v, sizeof(float)*chunk_size));
    CudaSafeCall(cudaMalloc( (void**)&data.w, sizeof(float)*chunk_size));
    CudaSafeCall(cudaMalloc( (void**)&data.data_real, sizeof(float)*chunk_size*nchan*nstokes));
    CudaSafeCall(cudaMalloc( (void**)&data.data_imag, sizeof(float)*chunk_size*nchan*nstokes));
    CudaSafeCall(cudaMalloc( (void**)&data.data_weight, sizeof(float)*chunk_size*nstokes));
    CudaSafeCall(cudaMalloc( (void**)&data.spw, sizeof(int)*chunk_size));
    CudaSafeCall(cudaMalloc( (void**)&data.field, sizeof(int)*chunk_size));
	size_t pb_size = sizeof(float)*nmaxcoords*nchan*nspw;
	cudaError err;
	err = cudaMalloc( (void**)&dev_coords.pb, pb_size);
	if(err != cudaSuccess)
	{
		if(err == cudaErrorMemoryAllocation)
		{
			fprintf(stderr, "Insuficient memory for primary beam data on device!");
			fprintf(stderr, "n_coords: %zu, nchan: %zu, nspw: %zu, size: %zu",
					nmaxcoords, nchan, nspw, pb_size);
		}
		else
		{
			fprintf(stderr, "Unknown error in allocation of dev_coords.pb (not cudaErrorMemoryAllocation)");
		}

		exit(-1);
	}
}/*}}}*/
void setup_freq(DataContainer& data, float* freq, const size_t nchan,/*{{{*/
                const size_t nspw)
{
    CudaSafeCall(cudaMalloc( (void**)&data.freq, sizeof(float)*nchan*nspw));
    CudaSafeCall(cudaMemcpy(data.freq, freq, sizeof(float)*nchan*nspw,
                cudaMemcpyHostToDevice));
}/*}}}*/
void copy_coords_to_cuda(/*{{{*/
                         Coords& coords, CoordContainer& dev_coords, 
                         float* freq, PrimaryBeam& pb, 
                         const int field, const size_t nchan,
                         const size_t nspw)
{
    dev_coords.n_coords = (size_t)coords.nStackPoints[field];

    float *omega = new float[dev_coords.n_coords*3];
    float *weight = new float[dev_coords.n_coords];
    float *pb_array = new float[dev_coords.n_coords*nchan*nspw];

    for(size_t coordID = 0; coordID < dev_coords.n_coords; coordID++)
    {
        omega[coordID*3+0] = coords.omega_x[field][coordID];
        omega[coordID*3+1] = coords.omega_y[field][coordID];
        omega[coordID*3+2] = coords.omega_z[field][coordID];
        weight[coordID]    = coords.weight[field][coordID];
        
        for(size_t spwID = 0; spwID < nspw; spwID++)
        {
            for(size_t chanID = 0; chanID < nchan; chanID++)
            {
                size_t index = spwID*nchan*dev_coords.n_coords
                             + chanID*dev_coords.n_coords
                             + coordID;
                pb_array[index] = pb.calc(coords.dx[field][coordID],
                                          coords.dy[field][coordID],
                                          freq[spwID*nchan+chanID]);
            }
        }
    }


	cudaError err;
	err = cudaMemcpyToSymbol( dev_omega, omega, 
			sizeof(float)*dev_coords.n_coords*3, 0,
			cudaMemcpyHostToDevice);
	if(err != cudaSuccess)
	{
		fprintf(stderr, "Failed to copy omega to __constant__ dev_omega");
		fprintf(stderr, "n_coords: %zu, size: %zu", dev_coords.n_coords,
				sizeof(float)*dev_coords.n_coords*3);
		exit(-1);
	}

    err = cudaMemcpyToSymbol( dev_weight, weight, 
                sizeof(float)*dev_coords.n_coords, 0,
                cudaMemcpyHostToDevice);
	if(err != cudaSuccess)
	{
		fprintf(stderr, "Failed to copy weight to __constant__ dev_weight");
		fprintf(stderr, "n_coords: %zu, size: %zu", dev_coords.n_coords,
				sizeof(float)*dev_coords.n_coords);
		exit(-1);
	}
	size_t pb_size = sizeof(float)*dev_coords.n_coords*nchan*nspw;
    err = cudaMemcpy( dev_coords.pb, pb_array, 
	                  pb_size, cudaMemcpyHostToDevice);
	if(err != cudaSuccess)
	{
		fprintf(stderr, "Failed to copy pb to device.");
		fprintf(stderr, "n_coords: %zu, size: %zu", dev_coords.n_coords,
				sizeof(float)*dev_coords.n_coords*nchan*nspw);
		exit(-1);
	}

    delete[] omega;
    delete[] weight;
    delete[] pb_array;
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
void visStack(DataContainer data, CoordContainer coords,/*{{{*/
              size_t chunk_size, size_t nchan, size_t nstokes)
{
	cuVisStack<<<BLOCKS,THREADS>>>(data, coords, chunk_size, nchan, nstokes);
}/*}}}*/

