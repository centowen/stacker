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

#ifdef CASACORE_VERSION_2
#include <casacore/casa/complex.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/VectorIter.h>
#include <casacore/casa/Arrays/MatrixIter.h>
#include <casacore/ms/MeasurementSets/MSTable.h>
#include <casacore/ms/MeasurementSets/MSColumns.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#else
#include <casa/complex.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/VectorIter.h>
#include <casa/Arrays/MatrixIter.h>
#include <ms/MeasurementSets/MSTable.h>
#include <ms/MeasurementSets/MSColumns.h>
#include <ms/MeasurementSets/MeasurementSet.h>
#endif

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

using casa::ROScalarColumn;

class Chunk;

class msio : public DataIO
{
	private:
		MatrixIterator<Complex>* datainit;
		VectorIterator<double>* uvwinit;
		VectorIterator<float>* weightinit;
		float* freq;
		size_t nchan, nspw;
		size_t nstokes;
		MeasurementSet* msin;
		MeasurementSet* msin_nonsorted;
		ROMSColumns* msincols;
		MeasurementSet* msout;
		MeasurementSet* msout_nonsorted;
		MSColumns* msoutcols;
		size_t currentVisibility;
		int readChunkDummy(Chunk& chunk);
		int readChunkSimple(Chunk& chunk);
		int readChunkIteratorbased(Chunk& chunk);

		int nfields;
		float* x_phase_centre;
		float* y_phase_centre;
		int datacolumn_;
		bool one_ptg_per_chunk_;
		int ptg_breaks_in_a_row;
		bool ptg_warning_done;

	public:
		static const int col_data = 0;
		static const int col_corrected_data = 1;
		static const int col_model_data = 2;

	public:
		msio(const char* msinfile, const char* msoutfile, 
		     int datacolumn = col_corrected_data,
			 const bool select_field = false, const char* field = "",
			 bool one_ptg_per_chunk = true);
		~msio();
		size_t nvis();

		int readChunk(Chunk& chunk);
		void writeChunk(Chunk& chunk);

		int nPointings();
		float xPhaseCentre(int fieldID);
		float yPhaseCentre(int fieldID);
		void setPhaseCentre(int fieldID, double x, double y);

		size_t nStokes();
		size_t nChan();
		size_t nSpw();
		float* getFreq(int spw);
};

#endif //inclusion guard
