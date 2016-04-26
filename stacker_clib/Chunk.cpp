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
	data_flag = NULL;
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
	data_flag_in = NULL;
	data_flag_out = NULL;
	weight_in = NULL;
	weight_out = NULL;

    update_datalinks();
}

Chunk::Chunk(const Chunk& c)
{
	dataset_id = c.dataset_id;
	nvis = c.nvis;
	max_nvis = c.nvis;

	inVis = new Visibility[this->nvis];
	outVis = new Visibility[this->nvis];

	nchan = c.nchan;
	nstokes = c.nstokes;

	if(nchan > 0 and nstokes > 0 and nvis > 0)
	{
		for(int i = 0; i < nvis; i++)
		{
			inVis[i].freq = c.inVis[i].freq;
			inVis[i].u = c.inVis[i].u;
			inVis[i].v = c.inVis[i].v;
			inVis[i].w = c.inVis[i].w;

			inVis[i].nstokes = c.inVis[i].nstokes;
			inVis[i].nchan = c.inVis[i].nchan;
			inVis[i].fieldID = c.inVis[i].fieldID;
			inVis[i].index = c.inVis[i].index;
			inVis[i].spw = c.inVis[i].spw;
			inVis[i].freq = c.inVis[i].freq;

			outVis[i].u = c.inVis[i].u;
			outVis[i].v = c.inVis[i].v;
			outVis[i].w = c.inVis[i].w;

			outVis[i].nstokes = c.inVis[i].nstokes;
			outVis[i].nchan = c.inVis[i].nchan;
			outVis[i].fieldID = c.inVis[i].fieldID;
			outVis[i].index = c.inVis[i].index;
			outVis[i].spw = c.inVis[i].spw;
		}

		data_real_in  = new float[nvis*nchan*nstokes];
		data_real_out = new float[nvis*nchan*nstokes];
		data_imag_in  = new float[nvis*nchan*nstokes];
		data_imag_out = new float[nvis*nchan*nstokes];
		data_flag_in  = new int[nvis*nchan*nstokes];
		data_flag_out = new int[nvis*nchan*nstokes];
		weight_in     = new float[nvis*nchan*nstokes];
		weight_out    = new float[nvis*nchan*nstokes];
		for(int i = 0; i < nchan*nstokes*nvis; i++)
		{
			data_real_in [i] = c.data_real_in [i];
			data_real_out[i] = c.data_real_out[i];
			data_imag_in [i] = c.data_imag_in [i];
			data_imag_out[i] = c.data_imag_out[i];
			data_flag_in [i] = c.data_flag_in [i];
			data_flag_out[i] = c.data_flag_out[i];
			weight_in    [i] = c.weight_in    [i];
			weight_out   [i] = c.weight_out   [i];
		}
	}

    update_datalinks();
}

Chunk::~Chunk()
{
	delete[] inVis;
	delete[] outVis;

	delete[] data_real_in;
	delete[] data_real_out;
	delete[] data_imag_in;
	delete[] data_imag_out;
	delete[] data_flag_in;
	delete[] data_flag_out;
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

void Chunk::resetSize()
{
	nvis = max_nvis;
    update_datalinks();
}

void Chunk::setSize(size_t size)
{
	if(size < 0)
		nvis = 0;
	else if(size <= max_nvis)
	{
		nvis = size;
	}
	else
	{
		std::cout << "setSize called with larger size!" << std::endl;
		delete[] inVis;
		delete[] outVis;
		inVis = new Visibility[size];
		outVis = new Visibility[size];
		nvis = size;
		max_nvis = size;
	}

    update_datalinks();
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
	delete[] data_flag_in;
	delete[] data_flag_out;
	delete[] weight_in;
	delete[] weight_out;
	data_real_in  = NULL;
	data_real_out = NULL;
	data_imag_in  = NULL;
	data_imag_out = NULL;
	data_flag_in  = NULL;
	data_flag_out = NULL;
	weight_in     = NULL;
	weight_out    = NULL;

	if(nchan > 0 and nstokes > 0 and nvis > 0)
	{
		data_real_in  = new float[max_nvis*nchan*nstokes];
		data_real_out = new float[max_nvis*nchan*nstokes];
		data_imag_in  = new float[max_nvis*nchan*nstokes];
		data_imag_out = new float[max_nvis*nchan*nstokes];
		data_flag_in  = new int[max_nvis*nchan*nstokes];
		data_flag_out = new int[max_nvis*nchan*nstokes];
		weight_in     = new float[max_nvis*nchan*nstokes];
		weight_out    = new float[max_nvis*nchan*nstokes];
		this->nchan  = nchan;
		this->nstokes = nstokes;
	}

    update_datalinks();
}

int Chunk::get_dataset_id()
{
	return dataset_id;
}

void Chunk::set_dataset_id(int id)
{
	dataset_id = id;
}

void Chunk::update_datalinks()
{
    if(nchan == (size_t)0 or nstokes == (size_t)0)
    {
        for(size_t i = 0; i < max_nvis; i++)
        {
            inVis[i].data_real  = NULL;
            inVis[i].data_imag  = NULL;
            inVis[i].data_flag  = NULL;
            inVis[i].weight     = NULL;
            outVis[i].data_real = NULL;
            outVis[i].data_imag = NULL;
            outVis[i].data_flag = NULL;
            outVis[i].weight    = NULL;
        }
    }

    if(nvis <= 0)
        return;

    for(size_t i = 0; i < nvis; i++)
    {
        inVis[i].data_real = &data_real_in[i*nchan*nstokes];
        inVis[i].data_imag = &data_imag_in[i*nchan*nstokes];
        inVis[i].data_flag = &data_flag_in[i*nchan*nstokes];
        inVis[i].weight    = &weight_in[i*nchan*nstokes];

        outVis[i].data_real = &data_real_out[i*nchan*nstokes];
        outVis[i].data_imag = &data_imag_out[i*nchan*nstokes];
        outVis[i].data_flag = &data_flag_out[i*nchan*nstokes];
        outVis[i].weight    = &weight_out[i*nchan*nstokes];
    }
}
