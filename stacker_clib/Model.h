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
#include "definitions.h"
#include "DataIO.h"

#ifndef __MODEL_H__
#define __MODEL_H__

using std::vector;
using std::ifstream;
using std::stringstream;
using std::sin;
using std::asin;
using std::cos;
using std::exp;
using std::log;
using std::cerr;
using std::cout;
using std::endl;
const int mod_point = 0;
const int mod_gaussian = 1;
const int mod_disk = 2;

class Model
{
private:
	struct stat statbuffer;
	string clfile;
	bool subtract_;
public:
	Model(string file, bool subtract);
	~Model();
	void compute(DataIO* ms, PrimaryBeam* pb);

public:
	int nPointings;
	int* nStackPoints;
	float** omega_x;
	float** omega_y;
	float** omega_z;

	int** model_type;
	float** omega_size;

	float** dx;
	float** dy;
	float** x;
	float** y;
	float** flux;
	float** size;
};

#endif

