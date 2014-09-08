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
#include "PrimaryBeam.h"
#include <iostream>
#include "definitions.h"
#include <limits.h>


PrimaryBeam::PrimaryBeam() {}
PrimaryBeam::~PrimaryBeam() {}
ConstantPrimaryBeam::ConstantPrimaryBeam() {}
ConstantPrimaryBeam::~ConstantPrimaryBeam() {}

ImagePrimaryBeam::ImagePrimaryBeam(const char fileName[])
{
	nx = 0;
	ny = 0;
	x0 = 0.;
	y0 = 0.;
	dx = 1.;
	dy = 1.;
	px_x0 = 0;
	px_y0 = 0;
	freq0 = 0.;
	data = NULL;
}

ImagePrimaryBeam::~ImagePrimaryBeam()
{
}


float ConstantPrimaryBeam::calc(float x, float y, float freq)
{
	return 1.;
}

// This function could be improved by adding interpolation,
// currently closest value is always used
float ImagePrimaryBeam::calc(float x, float y, float freq)
{
	float freqcomp = 1.;
	if(freq > tol && freq0 > tol)
		freqcomp = freq0/freq;

	int px_x, px_y;
	px_x = int(freqcomp*(x-x0)/dx+px_x0);
	px_y = int(freqcomp*(y-y0)/dy+px_y0);


	if(px_x < 0 || px_x >= nx) return 0.;
	if(px_y < 0 || px_y >= ny) return 0.;
	if(data == NULL) return 0.;
	return data[px_y][px_x];
}
