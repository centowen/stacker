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
#include "Chunk.h"
#include <iostream>

Visibility::Visibility()
{
	data_real = NULL;
	data_imag = NULL;
	weight = NULL;
	nstokes = 0;
	nchan = 0;
}

Visibility::~Visibility() {}

Chunk::Chunk(size_t size)
{
	dataset_id = dataset_none;
	nvis = size;
	max_nvis = size;

	inVis = new Visibility[size];
	outVis = new Visibility[size];

	nchan = 0;
	nstokes = 0;
	data_real_in = NULL;
	data_real_out = NULL;
	data_imag_in = NULL;
	data_imag_out = NULL;
	weight_in = NULL;
	weight_out = NULL;
}

Chunk::~Chunk()
{
	delete[] inVis;
	delete[] outVis;

	delete[] data_real_in;
	delete[] data_real_out;
	delete[] data_imag_in;
	delete[] data_imag_out;
	delete[] weight_in;
	delete[] weight_out;
	nvis = 0;
	nchan = 0;
	nstokes = 0;
}

size_t Chunk::size()
{
	return nvis;
}

void Chunk::setSize(int size)
{
	if(size < 0)
		nvis = 0;
	else if(size <= max_nvis)
	{
		nvis = size;
	}
	else
	{
		delete[] inVis, outVis;
		inVis = new Visibility[size];
		outVis = new Visibility[size];
		nvis = size;
		max_nvis = size;
	}
}

size_t Chunk::nChan()
{
	return nchan;
}

size_t Chunk::nStokes()
{
	return nstokes;
}

void Chunk::reshape_data(size_t nchan, size_t nstokes)
{
	// Calling this function with the same shape should be no-op.
	if(nchan == this->nchan and nstokes == this->nstokes)
		return;


	delete[] data_real_in;
	delete[] data_real_out;
	delete[] data_imag_in;
	delete[] data_imag_out;
	delete[] weight_in;
	delete[] weight_out;

	if(nchan > 0 and nstokes > 0 and nvis > 0)
	{
		data_real_in  = new float[nvis*nchan*nstokes];
		data_real_out = new float[nvis*nchan*nstokes];
		data_imag_in  = new float[nvis*nchan*nstokes];
		data_imag_out = new float[nvis*nchan*nstokes];
		weight_in     = new float[nvis*nchan*nstokes];
		weight_out    = new float[nvis*nchan*nstokes];
		this->nchan  = nchan;
		this->nstokes = nstokes;
	}

	for( int i = 0; i < nvis; i++)
	{
		inVis[i].data_real = &data_real_in[nchan*nstokes*i];
		inVis[i].data_imag = &data_imag_in[nchan*nstokes*i];
		inVis[i].weight    = &weight_in[nchan*nstokes*i];

		outVis[i].data_real = &data_real_out[nchan*nstokes*i];
		outVis[i].data_imag = &data_imag_out[nchan*nstokes*i];
		outVis[i].weight    = &weight_out[nchan*nstokes*i];
	}
}

int Chunk::get_dataset_id()
{
	return dataset_id;
}

void Chunk::set_dataset_id(int id)
{
	dataset_id = id;
}
