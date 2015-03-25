#ifndef __MODSUB_CHUNK_COMPUTER_GPU_CUDA_H__
#define __MODSUB_CHUNK_COMPUTER_GPU_CUDA_H__
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
#include "PrimaryBeam.h"
#include "Model.h"
#include "Chunk.h"
#include "DataIO.h"
#include "CommonCuda.h"

typedef struct _ModelContainer
{
    size_t n_mod_comp;
    float* pb; 
} ModelContainer;

void allocate_cuda_data_modsub(DataContainer& data, ModelContainer& dev_model,
                              const size_t nchan, const size_t nstokes,
                              const size_t chunk_size, const size_t nmax_mod_comp,
                              const size_t nspw);
void copy_model_to_cuda(Model& model, ModelContainer& dev_model, 
                         float* freq, PrimaryBeam& pb, 
                         const int field, const size_t nchan,
                         const size_t nspw);
void copy_model_to_cuda_partial(Model& model, ModelContainer& dev_model, 
                         float* freq, PrimaryBeam& pb, 
                         const int field, const size_t nchan,
                         const size_t nspw, 
						 const size_t n_mod_comp, const size_t first_mod_comp);
void modsub_chunk(DataContainer data, ModelContainer model, 
                  size_t chunk_size, size_t nchan, size_t n_stokes);

#endif // inclusion guard

