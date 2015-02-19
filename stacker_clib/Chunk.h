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

// #include <casa/complex.h>
// #include <casa/Arrays/Matrix.h>
// #include <casa/Arrays/Vector.h>
#include <iostream>

#ifndef __CHUNK_H__
#define __CHUNK_H__

// using casa::Complex;
// using casa::Matrix;
// using casa::Vector;

struct Visibility
{
	float u,v,w;
	float* freq;
	float *data_real, *data_imag, *weight;
	int nstokes, nchan;
	int fieldID, index, spw;

public:
	Visibility();
	~Visibility();
};

class Chunk
{
public:
	static const int dataset_none = -1;


private:
	int dataset_id;
	size_t nvis, max_nvis;
	size_t nchan;
	size_t nstokes;

public:
	float* data_real_in;
	float* data_imag_in;
	float* weight_in;
	float* data_real_out;
	float* data_imag_out;
	float* weight_out;
	Visibility *inVis, *outVis;

	Chunk(size_t size);
	~Chunk();

	size_t size();
	void setSize(int size);

	size_t nChan();
	size_t nStokes();
	void reshape_data(size_t nchan, size_t nstokes);

	int get_dataset_id();
	void set_dataset_id(int id);
};

#endif // end of inclusion guard
