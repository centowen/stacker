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

#ifndef __DEFINITIONS_H__
#define __DEFINITIONS_H__

const double tol = 1e-9;
const int N_THREAD = 24;
const int N_CHUNK = 2*N_THREAD;
const int CHUNK_SIZE = 10000;
const double c = 299792458;

const int PB_CONST = 0;
const int PB_MS = 1;
const int PB_FITS = 2;

const int FILE_TYPE_NONE = 0;
const int FILE_TYPE_MS = 1;
const int FILE_TYPE_FITS = 2;

// Operate on the data column rather than the corrected data column.
const int MS_DATACOLUMN_DATA = 1;

#endif
