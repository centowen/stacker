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

Visibility::~Visibility()
{
	if(nstokes>0 && nchan > 0)
	{
		delete[] data_real, data_imag, weight;
	}
}

Chunk::Chunk(int size)
{
	inVis = new Visibility[size];
	outVis = new Visibility[size];
	nvis = size;
}
Chunk::~Chunk()
{
	delete[] inVis;
	delete[] outVis;
	nvis = 0;
}

int Chunk::size()
{
	return nvis;
}

void Chunk::setSize(int size)
{
	if(size < 0)
		nvis = 0;
	else
	{
// 		delete[] inVis, outVis;
// 		inVis = new Visibility[size];
// 		outVis = new Visibility[size];
		nvis = size;
	}
}
