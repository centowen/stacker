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

#include "MSComputer.h"
#include "Coords.h"
#include "Model.h"
#include "PrimaryBeam.h"
#include "DataIO.h"
#include <pthread.h>

#include "StackChunkComputerGpu_cuda.h"
#include "ModsubChunkComputerGpu_cuda.h"
#include "StackMCCCGpu_cuda.h"

#ifndef __STACK_MC_CC_GPU_H__
#define __STACK_MC_CC_GPU_H__

class StackMCCCGpu: public ChunkComputer
{
	private:
		Coords** coords;
		Model** models;
		PrimaryBeam* pb;

		DataContainer dev_data;
		ModelContainer dev_model;
		CoordContainer dev_coords;
		MCResultContainer dev_results;

		int nbin;
		int nmc;
		float* bins;
		float* res_flux;
		float* res_weight;

		int field;
		float* freq;
		int nspw;

		bool redoWeights;

	public:
		StackMCCCGpu(Coords** coords, Model** models, PrimaryBeam* pb, 
		             int nmc, double* bins, int nbin);

		// Called from computer and allows to access data,
		// unlike normal constructor which is called before computer
		// is created.
		void preCompute(DataIO* ms);
		virtual void computeChunk(Chunk* chunk);
		void postCompute(DataIO* ms);

		double* get_flux();
		double* get_weight();
};

#endif // inclusion guard

