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
#include "CachedDataIO.h"
#include "definitions.h"

CachedDataIO::CachedDataIO(DataIO* dataio, size_t max_chunks) : DataIO(), max_chunks_(max_chunks)
{
	nfields = dataio->nPointings();
	x_phase_centre = new float[nfields];
	y_phase_centre = new float[nfields];

	for(int fieldID = 0; fieldID < nfields; fieldID++)
	{
		x_phase_centre[fieldID] = dataio->xPhaseCentre(fieldID);
		y_phase_centre[fieldID] = dataio->yPhaseCentre(fieldID);
	}

	n_vis   = dataio->nvis();
	nchan   = dataio->nChan();
	nspw    = dataio->nSpw();
	nstokes = dataio->nStokes();

	bool chunks_left = true;
	while(chunks_left and cache.size() < max_chunks_)
	{
		cache.push_back(Chunk(CHUNK_SIZE));
		if(dataio->readChunk(*cache.rbegin()) <= 0)
		{
			chunks_left = false;
			cache.pop_back();
		}
	}

	cache_iterator = cache.begin();

	freq = new float[nchan*nspw];
	for(size_t spw = 0; spw < nspw; spw++)
	{
		std::copy(dataio->getFreq(spw), &dataio->getFreq(spw)[nchan], freq);
	}
}

size_t CachedDataIO::nvis()
{
	return n_vis;
}

void CachedDataIO::writeChunk(Chunk& chunk)/*{{{*/
{ 
} /*}}}*/

size_t CachedDataIO::readChunk(Chunk& chunk)/*{{{*/
{
	// FIXME: I have no idea what I am doing here !!!
	if(cache_iterator == cache.end())
		return 0;

	chunk.set_dataset_id(dataset_id);
	chunk.reshape_data(this->nchan, this->nstokes);

	Chunk& chunk_internal = *cache_iterator;
	cache_iterator++;

	size_t len = chunk_internal.size()*chunk_internal.nChan()*chunk_internal.nStokes();

	// Copy all data from the old to the new chunk.
	// This copy should perhaps be done in Chunk.
	std::copy(chunk_internal.data_real_in,
			  &chunk_internal.data_real_in[len],
			  chunk.data_real_in);
	std::copy(chunk_internal.data_imag_in,
			  &chunk_internal.data_imag_in[len],
			  chunk.data_imag_in);
	std::copy(chunk_internal.weight_in,
			  &chunk_internal.weight_in[len],
			  chunk.weight_in);
	std::copy(chunk_internal.data_real_out,
			  &chunk_internal.data_real_out[len],
			  chunk.data_real_out);
	std::copy(chunk_internal.data_imag_out,
			  &chunk_internal.data_imag_out[len],
			  chunk.data_imag_out);
	std::copy(chunk_internal.weight_out,
			  &chunk_internal.weight_out[len],
			  chunk.weight_out);
	std::copy(chunk_internal.inVis,
			  &chunk_internal.inVis[chunk_internal.size()],
			  chunk.inVis);
	std::copy(chunk_internal.outVis,
			  &chunk_internal.outVis[chunk_internal.size()],
			  chunk.outVis);

	chunk.update_datalinks();

	return 0;
}/*}}}*/

int CachedDataIO::nPointings()
{
	return nfields;
}

float CachedDataIO::xPhaseCentre(int id)
{
	return x_phase_centre[id];
}

float CachedDataIO::yPhaseCentre(int id)
{
	return y_phase_centre[id];
}

void CachedDataIO::setPhaseCentre(int fieldID, double x, double y)
{
}

size_t CachedDataIO::nChan()
{
	return nchan;
}

size_t CachedDataIO::nSpw()
{
	return nspw;
}

size_t CachedDataIO::nStokes()
{
	return nstokes;
}

float* CachedDataIO::getFreq(int spw)
{
	return &freq[spw*nchan];
}

void CachedDataIO::restart()
{
	cache_iterator = cache.begin();
}
