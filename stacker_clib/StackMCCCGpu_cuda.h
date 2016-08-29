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
#ifndef __STACK_MC_CC_GPU_CUDA_H__
#define __STACK_MC_CC_GPU_CUDA_H__
typedef struct _MC_ResultContainer
{
	float* bins;
	float* flux;
	float* weight;
	int nbin;
} MCResultContainer;

void allocate_cuda_data_stack_mc(MCResultContainer& results, float* bins, int nbin);
void zero_results_stack_mc(MCResultContainer& results, float* bins);
void compute_results_stack_mc(MCResultContainer& results, DataContainer& data,
                                         size_t chunk_size, size_t nchan,
                                         size_t nstokes);
void copy_results_to_host(MCResultContainer& results, float* flux,
                          float* weight, int mc_id);
#endif // inclusion guard
