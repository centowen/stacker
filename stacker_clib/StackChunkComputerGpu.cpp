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

#include "StackChunkComputerGpu.h"
#include "Chunk.h"
#include "Coords.h"
#include "PrimaryBeam.h"

StackChunkComputerGpu::StackChunkComputerGpu(Coords* coords, PrimaryBeam* pb)/*{{{*/
{
	this->coords = coords;
	this->pb = pb;
}/*}}}*/
void StackChunkComputerGpu::computeChunk(Chunk* chunk) /*{{{*/
{
	copy_data_to_cuda(dev_data, *chunk);
	for(size_t uvrow = 0; uvrow < chunk->size(); uvrow++)
	{
		Visibility& inVis = chunk->inVis[uvrow];
		Visibility& outVis = chunk->outVis[uvrow];
		if(coords->nStackPoints[inVis.fieldID] > 0)
			outVis.fieldID = 0;
		else
			outVis.fieldID = 1;
	}
	visStack(dev_data, dev_coords, chunk->size(),
	         chunk->nChan(), chunk->nStokes());
	copy_data_to_host(dev_data, *chunk);
}/*}}}*/
void StackChunkComputerGpu::preCompute(DataIO* dataio)/*{{{*/
{
	// If we only have one field in the data the weights do not need to be 
	// updated. 
	if(dataio->nPointings() > 1)
	{
		redoWeights = true;
	}


	pthread_mutex_init(&fluxMutex, NULL);
    sumvisweight = 0.;
    sumweight = 0.;

	coords->computeCoords(dataio, *pb);

	allocate_cuda_data(dev_data, dataio->nChan(),
	                   dataio->nStokes(), CHUNK_SIZE);
	setup_freq(dev_data, *dataio);
	copy_coords_to_cuda(*coords, dev_coords, *dataio, *pb);
}/*}}}*/
void StackChunkComputerGpu::postCompute(DataIO* data)/*{{{*/
{
	// After stacking we have all visibilities in field 0.
	// Set centre to 0., 0. to indicate that coordinates at this point are
	// arbitrary.
	data->setPhaseCentre(0, 0., 0.);

	// Maybe it would be nice to also remove all the other fields?
	// Also should properly flag visibilities that got bad in stacking.
}/*}}}*/
double StackChunkComputerGpu::flux()/*{{{*/
{
    return 0.;
}/*}}}*/
