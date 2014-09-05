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

#include <pthread.h>
#include <casa/complex.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/VectorIter.h>
#include <casa/Arrays/MatrixIter.h>
#include <ms/MeasurementSets/MSTable.h>
#include <ms/MeasurementSets/MSColumns.h>
#include <ms/MeasurementSets/MeasurementSet.h>

#include "DataIO.h"

#ifndef __MSIO_H__
#define __MSIO_H__

using casa::MeasurementSet;
using casa::MSColumns;
using casa::ROMSColumns;
using casa::MatrixIterator;
using casa::VectorIterator;
using casa::Array;
using casa::Complex;
using casa::Matrix;
using casa::Vector;
using casa::Float;
using casa::IPosition;

class Chunk;

class msio : public DataIO
{
	private:
		pthread_mutex_t mutex;
		MatrixIterator<Complex>* datainit;
		VectorIterator<double>* uvwinit;
		VectorIterator<float>* weightinit;
		double** freq;
		int nchan, nspw;
		MeasurementSet* msin;
		ROMSColumns* msincols;
		MeasurementSet* msout;
		MSColumns* msoutcols;
		int currentVisibility;
		int readChunkDummy(Chunk& chunk);
		int readChunkSimple(Chunk& chunk);
		int readChunkIteratorbased(Chunk& chunk);

		int nfields;
		float* x_phase_centre;
		float* y_phase_centre;
		int datacolumn;

		static const int col_data = 0;
		static const int col_corrected_data = 1;

	public:
		msio(const char* msinfile, const char* msoutfile, 
		          pthread_mutex_t* mutex, bool workondatacolumn = false);
		~msio();
		int nvis();

		int readChunk(Chunk& chunk);
		void writeChunk(Chunk& chunk);

		int nPointings();
		float xPhaseCentre(int fieldID);
		float yPhaseCentre(int fieldID);
		void setPhaseCentre(int fieldID, double x, double y);

};

#endif //inclusion guard
