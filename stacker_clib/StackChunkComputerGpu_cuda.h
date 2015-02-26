#ifndef __STACK_CHUNK_COMPUTER_GPU_CUDA_H__
#define __STACK_CHUNK_COMPUTER_GPU_CUDA_H__
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
#include "Coords.h"
#include "Chunk.h"
#include "DataIO.h"
// class PrimaryBeam;
// class Coords;
// class Chunk;
// class DataIO;

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

typedef struct _CoordContainer
{
    size_t n_coords;
    float* pb; 
} CoordContainer;

void allocate_cuda_data(DataContainer& data, CoordContainer& dev_coords,
                        const size_t nchan, const size_t nstokes,
                        const size_t chunk_size, const size_t nmaxcoords,
                        const size_t nspw);
void setup_freq(DataContainer& data, float* freq,
                const size_t nchan, const size_t nspw);
void copy_coords_to_cuda(Coords& coords, CoordContainer& dev_coords, 
                         float* freq, PrimaryBeam& pb, 
                         const int field, const size_t nchan,
                         const size_t nspw);
void copy_data_to_cuda(DataContainer& data, Chunk& chunk);
void copy_data_to_host(DataContainer& data, Chunk& chunk);
void visStack(DataContainer data, CoordContainer coords, 
              size_t chunk_size, size_t nchan, size_t n_stokes);

#endif // inclusion guard

