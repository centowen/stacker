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

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <cmath>

#include <sys/stat.h>

#include "PrimaryBeam.h"
#include "DataIO.h"
#include "definitions.h"

#ifndef __COORDS_H__
#define __COORDS_H__

using std::vector;
using std::ifstream;
using std::stringstream;
using std::sin;
using std::asin;
using std::cos;
using std::cerr;
using std::cout;
using std::endl;

#define pi M_PI

class Coords
{
private:
	struct stat statbuffer;
public:
	Coords(const char* coordfile);
	Coords(double* x, double* y, double* weight, int nstack);
	~Coords();
	void computeCoords(DataIO* ms, PrimaryBeam& pb);

public:
	int nPointings;
	float** omega_x;
	float** omega_y;
	float** omega_z;
	float** dx;
	float** dy;
	double* x_raw;
	double* y_raw;
	double* weight_raw;
	int nStackPoints_raw;
	float** x;
	float** y;
	float** weight;
	int* nStackPoints;
	int nStackPointsVisible;
};

#endif
