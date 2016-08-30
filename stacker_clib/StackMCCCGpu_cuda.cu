// stacker, Python module for stacking of interferometric data.
// Copyright (C) 2014  Lukas Lindroos
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// // This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. 
//
// Library to stack and modsub ms data.
#include <cuda.h>
#include "CommonCuda.h"
#include "StackMCCCGpu_cuda.h"
#include "Chunk.h"
#include "definitions.h"
#include "cuda_error.h"
#include <iostream>
using std::cout;
using std::endl;

void allocate_cuda_data_stack_mc(MCResultContainer& results, float* bins, int nbin)/*{{{*/
{
	CudaSafeCall(cudaMalloc( (void**)&results.bins, sizeof(float)*(nbin+1)));
	CudaSafeCall(cudaMalloc( (void**)&results.flux, sizeof(float)*(nbin)));
	CudaSafeCall(cudaMalloc( (void**)&results.weight, sizeof(float)*(nbin)));

	float* zeros = new float[nbin];
	for(int i = 0; i < nbin; i++)
		zeros[i] = 0.;

    CudaSafeCall(cudaMemcpy(results.bins, bins, sizeof(float)*(nbin+1),
                            cudaMemcpyHostToDevice));
    CudaSafeCall(cudaMemcpy(results.flux, zeros, sizeof(float)*(nbin),
                            cudaMemcpyHostToDevice));
    CudaSafeCall(cudaMemcpy(results.weight, zeros, sizeof(float)*(nbin),
                            cudaMemcpyHostToDevice));
	delete[] zeros;
	results.nbin = nbin;
};/*}}}*/
void zero_results_stack_mc(MCResultContainer& results, float* bins)/*{{{*/
{
	float* zeros = new float[results.nbin];
	for(int i = 0; i < results.nbin; i++)
		zeros[i] = 0.;

    CudaSafeCall(cudaMemcpy(results.bins, bins, sizeof(float)*(results.nbin+1),
                            cudaMemcpyHostToDevice));
    CudaSafeCall(cudaMemcpy(results.flux, zeros, sizeof(float)*(results.nbin),
                            cudaMemcpyHostToDevice));
    CudaSafeCall(cudaMemcpy(results.weight, zeros, sizeof(float)*(results.nbin),
                            cudaMemcpyHostToDevice));
	delete[] zeros;
}/*}}}*/
__global__ void cu_compute_results_stack_mc(MCResultContainer results, DataContainer data,/*{{{*/
		                                 size_t chunk_size, size_t nchan, size_t nstokes)
{
    size_t uvrow = threadIdx.x + blockIdx.x*blockDim.x;
    while(uvrow < chunk_size)
    {
		int uvbin = 0;
		float uvdist = sqrt(data.u[uvrow]*data.u[uvrow] + data.v[uvrow]*data.v[uvrow]);
		for(size_t i = 0; i < results.nbin+1; i++)
		{
			if(uvdist < results.bins[i])
			{
				uvbin = i;
				break;
			}
		}
// 		while((uvbin <= results.nbin) && (uvdist < results.bins[uvbin]))
// 			uvbin++;
// 		while( (uvbin <= results.nbin) && (uvdist2 < (results.bins[uvbin]*results.bins[uvbin])))
// 			uvbin++;
		if(uvbin > 0 && uvbin <= results.nbin)
		{
			uvbin -= 1;
			float weight = 0.;
			float weighteddata = 0.;
			for(size_t chanID = 0; chanID < nchan; chanID++)
			{
				size_t stokesID = 0;
				size_t weightindex = uvrow*nstokes+stokesID;
				size_t dataindex = nchan*weightindex + chanID;
				if(!data.data_flag[dataindex])
				{
					weighteddata += data.data_real[dataindex]*
									data.data_weight[weightindex];
					weight += data.data_weight[weightindex];
				}
// 				weighteddata += data.data_real[dataindex]*
// 								data.data_weight[weightindex];
// 				weight += data.data_weight[weightindex];
			}
			atomicAdd(&results.flux[uvbin],
					  weighteddata);
			atomicAdd(&results.weight[uvbin],
					  weight);
		}
        uvrow += blockDim.x*gridDim.x;
	}
};/*}}}*/
void compute_results_stack_mc(MCResultContainer& results, /*{{{*/
                                         DataContainer& data, size_t chunk_size,
                                         size_t nchan, size_t nstokes)
{
	cu_compute_results_stack_mc<<<BLOCKS,THREADS>>>(results, data, chunk_size, 
	                                                nchan, nstokes);
}/*}}}*/
void copy_results_to_host(MCResultContainer& results, float* flux,/*{{{*/
		                  float* weight, int mc_id)
{
	float* flux_buff = new float[results.nbin];
	float* weight_buff = new float[results.nbin];

    CudaSafeCall(cudaMemcpy(flux_buff, results.flux,
                sizeof(float)*results.nbin,
                cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaMemcpy(weight_buff, results.weight,
                sizeof(float)*results.nbin,
                cudaMemcpyDeviceToHost));

	for(int bin = 0; bin < results.nbin; bin++)
	{
// 		flux[mc_id*results.nbin+bin] = 1.;
// 		weight[mc_id*results.nbin+bin] = 1.;
		flux[mc_id*results.nbin+bin] += flux_buff[bin];
		weight[mc_id*results.nbin+bin] += weight_buff[bin];
	}

	delete[] flux_buff;
	delete[] weight_buff;
};/*}}}*/
