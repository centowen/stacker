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
#include <iostream>
#include "ModsubChunkComputerGpu_cuda.h"
#include "DataIO.h"
#include "PrimaryBeam.h"
#include "Chunk.h"
#include "Model.h"
#include "definitions.h"
#include "cuda_error.h"

__constant__ float dev_omega[3*N_MAX_MOD_COMP];
__constant__ float dev_omega_size[N_MAX_MOD_COMP];
__constant__ float dev_flux[N_MAX_MOD_COMP];

__global__ void cu_modsub(DataContainer data, ModelContainer model, /*{{{*/
                          size_t chunk_size, size_t nchan, size_t nstokes)
{
    size_t uvrow = threadIdx.x + blockIdx.x*blockDim.x;

    float exponent, sin_exponent, cos_exponent, pbcor;
	float m_real, m_imag;
	float extent;

    while(uvrow < chunk_size)
    {
        float* freq = &data.freq[data.spw[uvrow]*nchan];

        for(size_t chanID = 0; chanID < nchan; chanID++)
        {
            m_real = 0.;
            m_imag = 0.;

            for(size_t compID = 0; compID < model.n_mod_comp; compID++)
            {
                size_t pbindex = data.spw[uvrow]*nchan*model.n_mod_comp
                               + chanID*model.n_mod_comp
                               + compID;
                pbcor = model.pb[pbindex];
                exponent = freq[chanID]*(dev_omega[compID*3+0]*data.u[uvrow] + 
                                         dev_omega[compID*3+1]*data.v[uvrow] + 
                                         dev_omega[compID*3+2]*data.w[uvrow]);
                sincos(exponent, &sin_exponent, &cos_exponent);
// 				cos_exponent = 1.;
// 				sin_exponent = 0.;
// 				extent = 1.;
				if(dev_omega_size[compID] > 1e-10)
				{
					extent = exp(-freq[chanID]*freq[chanID] *
					              (data.u[uvrow]*data.u[uvrow] +
					               data.v[uvrow]*data.v[uvrow]) *
					              dev_omega_size[compID]);
				}
				else
				{
					extent = 1.;
				}
                m_real += pbcor*dev_flux[compID]*extent*cos_exponent;
                m_imag += pbcor*dev_flux[compID]*extent*sin_exponent;
            }

            for(size_t stokesID = 0; stokesID < nstokes; stokesID++)
            {
                size_t weightindex = uvrow*nstokes+stokesID;
                size_t dataindex = nchan*weightindex + chanID;

// 				data.data_real[dataindex] -= m_real;
// 				data.data_imag[dataindex] -= m_imag;
				data.data_real[dataindex] = data.data_real[dataindex] - m_real;
				data.data_imag[dataindex] = data.data_imag[dataindex] - m_imag;
            }
        }

        uvrow += blockDim.x*gridDim.x;
    }
};/*}}}*/

void allocate_cuda_data_modsub(DataContainer& data, ModelContainer& dev_model,/*{{{*/
                               const size_t nchan, const size_t nstokes,
                               const size_t chunk_size, const size_t nmax_mod_comp,
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
	size_t pb_size = sizeof(float)*nmax_mod_comp*nchan*nspw;
	cudaError err;
	err = cudaMalloc( (void**)&dev_model.pb, pb_size);
	if(err != cudaSuccess)
	{
		if(err == cudaErrorMemoryAllocation)
		{
			fprintf(stderr, "Insuficient memory for primary beam data on device!\n");
			fprintf(stderr, "n_mod_comp: %zu, nchan: %zu, nspw: %zu, size: %zu\n",
					nmax_mod_comp, nchan, nspw, pb_size);
		}
		else
		{
			fprintf(stderr, "Unknown error in allocation of dev_model.pb (not cudaErrorMemoryAllocation)\n");
			fprintf(stderr, "Error: %s\n", cudaGetErrorString(err));
		}

		exit(-1);
	}
}/*}}}*/
void copy_model_to_cuda(/*{{{*/
                        Model& model, ModelContainer& dev_model, 
                        float* freq, PrimaryBeam& pb, 
                        const int field, const size_t nchan,
                        const size_t nspw)
{
    dev_model.n_mod_comp = (size_t)model.nStackPoints[field];

    float *omega = new float[dev_model.n_mod_comp*3];
    float *omega_size = new float[dev_model.n_mod_comp];
    float *flux = new float[dev_model.n_mod_comp];
    float *pb_array = new float[dev_model.n_mod_comp*nchan*nspw];

    for(size_t mod_comp_id = 0; mod_comp_id < dev_model.n_mod_comp; mod_comp_id++)
    {
        omega[mod_comp_id*3+0]  = model.omega_x   [field][mod_comp_id];
        omega[mod_comp_id*3+1]  = model.omega_y   [field][mod_comp_id];
        omega[mod_comp_id*3+2]  = model.omega_z   [field][mod_comp_id];
        omega_size[mod_comp_id] = model.omega_size[field][mod_comp_id];
        flux[mod_comp_id]       = model.flux      [field][mod_comp_id];
        
        for(size_t spwID = 0; spwID < nspw; spwID++)
        {
            for(size_t chanID = 0; chanID < nchan; chanID++)
            {
                size_t index = spwID*nchan*dev_model.n_mod_comp
                             + chanID*dev_model.n_mod_comp
                             + mod_comp_id;
                pb_array[index] = pb.calc(model.dx[field][mod_comp_id],
                                          model.dy[field][mod_comp_id],
                                          freq[spwID*nchan+chanID]);
            }
        }
    }


	cudaError err;
	err = cudaMemcpyToSymbol( dev_omega, omega, 
			sizeof(float)*dev_model.n_mod_comp*3, 0,
			cudaMemcpyHostToDevice);
	if(err != cudaSuccess)
	{
		fprintf(stderr, "Failed to copy omega to __constant__ dev_omega\n");
		fprintf(stderr, "n_mod_comp: %zu, size: %zu\n", dev_model.n_mod_comp,
				sizeof(float)*dev_model.n_mod_comp*3);
		exit(-1);
	}

	err = cudaMemcpyToSymbol( dev_omega_size, omega_size, 
			sizeof(float)*dev_model.n_mod_comp, 0,
			cudaMemcpyHostToDevice);
	if(err != cudaSuccess)
	{
		fprintf(stderr, "Failed to copy omega_size to __constant__ dev_omega_size\n");
		fprintf(stderr, "n_mod_comp: %zu, size: %zu\n", dev_model.n_mod_comp,
				sizeof(float)*dev_model.n_mod_comp);
		exit(-1);
	}

    err = cudaMemcpyToSymbol( dev_flux, flux, 
                sizeof(float)*dev_model.n_mod_comp, 0,
                cudaMemcpyHostToDevice);
	if(err != cudaSuccess)
	{
		fprintf(stderr, "Failed to copy flux to __constant__ dev_flux\n");
		fprintf(stderr, "n_mod_comp: %zu, size: %zu\n", dev_model.n_mod_comp,
				sizeof(float)*dev_model.n_mod_comp);
		exit(-1);
	}
	size_t pb_size = sizeof(float)*dev_model.n_mod_comp*nchan*nspw;
    err = cudaMemcpy( dev_model.pb, pb_array, 
	                  pb_size, cudaMemcpyHostToDevice);
	if(err != cudaSuccess)
	{
		fprintf(stderr, "Failed to copy pb to device.\n");
		fprintf(stderr, "n_model: %zu, size: %zu\n", dev_model.n_mod_comp,
				sizeof(float)*dev_model.n_mod_comp*nchan*nspw);
		exit(-1);
	}

    delete[] omega;
    delete[] omega_size;
    delete[] flux;
    delete[] pb_array;
}/*}}}*/
void modsub_chunk(DataContainer data, ModelContainer model,/*{{{*/
              size_t chunk_size, size_t nchan, size_t nstokes)
{
	cu_modsub<<<BLOCKS,THREADS>>>(data, model, chunk_size, nchan, nstokes);
}/*}}}*/

