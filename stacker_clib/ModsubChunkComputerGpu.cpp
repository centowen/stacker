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

#include "ModsubChunkComputerGpu.h"
#include "Chunk.h"
#include "Model.h"
#include "PrimaryBeam.h"

ModsubChunkComputerGpu::ModsubChunkComputerGpu(Model* model, PrimaryBeam* pb)/*{{{*/
{
	this->model = model;
	this->pb = pb;
	field = -999;
	freq = NULL;
	nspw = 0;

}/*}}}*/

void ModsubChunkComputerGpu::computeChunk(Chunk* chunk) /*{{{*/
{
	if(chunk->inVis[0].fieldID != field)
	{
		field = chunk->inVis[0].fieldID;
		copy_model_to_cuda(*model, dev_model, freq, *pb, field, chunk->nChan(), nspw);
	}
	copy_data_to_cuda(dev_data, *chunk);
	for(size_t uvrow = 0; uvrow < chunk->size(); uvrow++)
	{
		Visibility& inVis = chunk->inVis[uvrow];
		Visibility& outVis = chunk->outVis[uvrow];
		outVis.fieldID = inVis.fieldID;
	}
	modsub_chunk(dev_data, dev_model, chunk->size(),
	         chunk->nChan(), chunk->nStokes());
	copy_data_to_host(dev_data, *chunk);
}/*}}}*/
void ModsubChunkComputerGpu::preCompute(DataIO* dataio)/*{{{*/
{
	model->compute(dataio, pb);

	int nmax_model_comp = 0;
	for(int i = 0; i < model->nPointings; i++)
		if(model->nStackPoints[i] > nmax_model_comp)
		{
			nmax_model_comp = model->nStackPoints[i];
		}

	allocate_cuda_data_modsub(dev_data, dev_model, dataio->nChan(),
	                          dataio->nStokes(), CHUNK_SIZE, 
					          nmax_model_comp, dataio->nSpw());

	nspw = dataio->nSpw();
    freq = new float[dataio->nChan()*dataio->nSpw()];
    // Load frequencies into freq[].
    for(size_t chanID = 0; chanID < dataio->nChan(); chanID++)
    {   
        for(size_t spwID = 0; spwID < dataio->nSpw(); spwID++)
        {   
            freq[spwID*dataio->nChan()+chanID] = (float)dataio->getFreq(spwID)[chanID];
        }   
    }   
	// Copy freq to device.
	setup_freq(dev_data, freq, dataio->nChan(), dataio->nSpw());
}/*}}}*/
void ModsubChunkComputerGpu::postCompute(DataIO* data)/*{{{*/
{
	// FIXME: Should clean up cuda here.
	delete[] freq;
}/*}}}*/
