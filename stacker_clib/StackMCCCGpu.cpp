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
#include <iostream>

#include "StackMCCCGpu.h"
#include "Chunk.h"
#include "Coords.h"
#include "Model.h"
#include "PrimaryBeam.h"


StackMCCCGpu::StackMCCCGpu(Coords** coordlists, Model** models, PrimaryBeam* pb,/*{{{*/
                           int nmc, double* bins, int nbin)
{
	this->coords = coordlists;
	this->models = models;
	this->pb = pb;
	field = -999;
	freq = NULL;
	nspw = 0;
	this->nmc = nmc;
	this->nbin = nbin;
	this->bins = new float[nbin];
	for(int i = 0; i < nbin+1; i++)
		this->bins[i] = bins[i];
	res_flux = new float[nmc*nbin];
	res_weight = new float[nmc*nbin];
	for(int i = 0; i < nbin*nmc; i++)
	{
		res_flux[i] = 0.;
		res_weight[i] = 0.;
	}
}/*}}}*/
void StackMCCCGpu::computeChunk(Chunk* chunk)/*{{{*/
{

	for(int i = 0; i < nmc; i++)
	{
// 		Chunk local_chunk(*chunk);
		field = chunk->inVis[0].fieldID;
		copy_coords_to_cuda(*coords[i], dev_coords, freq, *pb, field, chunk->nChan(), nspw);
		copy_data_to_cuda(dev_data, *chunk);
		visStack(dev_data, dev_coords, chunk->size(), 
				 chunk->nChan(), chunk->nStokes());
// 		copy_data_to_host(dev_data, local_chunk);
		zero_results_stack_mc(dev_results);
		compute_results_stack_mc(dev_results, dev_data, chunk->size(),
		                         chunk->nChan(), chunk->nStokes());
		copy_results_to_host(dev_results, res_flux, res_weight, i);
	}
}/*}}}*/
void StackMCCCGpu::preCompute(DataIO* dataio)/*{{{*/
{
	if(dataio->nPointings() > 1)
	{
		redoWeights = true;
	}

	int nmax_coords = 0;
	int nmax_model_comp = 0;

	for(int i = 0; i < nmc; i++)
	{
		coords[i]->computeCoords(dataio, *pb);
		for(int i_ptg = 0; i_ptg < coords[i]->nPointings; i_ptg++)
			if(coords[i]->nStackPoints[i_ptg] > nmax_coords)
				nmax_coords = coords[i]->nStackPoints[i_ptg];


		if(models[i] != NULL)
		{
			models[i]->compute(dataio, pb);
			for(int i_ptg = 0; i_ptg < models[i]->nPointings; i_ptg++)
				if(models[i]->nStackPoints[i_ptg] > nmax_model_comp)
					nmax_model_comp = models[i]->nStackPoints[i_ptg];
		}
	}

	allocate_cuda_data(dev_data, dataio->nChan(), 
	                   dataio->nStokes(), CHUNK_SIZE);
	allocate_cuda_data_stack(dev_coords, dataio->nChan(),
							 nmax_coords, dataio->nSpw());
	allocate_cuda_data_stack_mc(dev_results, bins, nbin);

	nspw = dataio->nSpw();
	freq = new float[dataio->nChan()*dataio->nSpw()];

	for(size_t chanID = 0; chanID < dataio->nChan(); chanID++)
	{
		for(size_t spwID = 0; spwID < nspw; spwID++)
		{
			freq[spwID*dataio->nChan()+chanID] = (float)dataio->getFreq(spwID)[chanID];
		}
	}
	setup_freq(dev_data, freq, dataio->nChan(), dataio->nSpw());

}/*}}}*/
void StackMCCCGpu::postCompute(DataIO* data)/*{{{*/
{
	delete[] freq;
}/*}}}*/
double* StackMCCCGpu::get_flux()/*{{{*/
{
	double* ret = new double[nmc*nbin];
	for(int i = 0; i < nmc*nbin; i++)
	{
		ret[i] = (double)res_flux[i];
	}
	return ret;
}/*}}}*/
double* StackMCCCGpu::get_weight()/*{{{*/
{
	double* ret = new double[nmc*nbin];
	for(int i = 0; i < nmc*nbin; i++)
	{
		ret[i] = (double)res_weight[i];
	}
	return ret;
}/*}}}*/
