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
#include <fitsio.h>
#include <exception>
#include <stdexcept>

using std::string;

#include "DataIO.h"

#ifndef __DATA_IO_FITSIO__
#define __DATA_IO_FITSIO__

class Chunk;

class keynameError : public std::runtime_error
{
	public:
		keynameError(const char* m = "") throw (): std::runtime_error(m) { }
		~keynameError() throw () {};
};

class DataIOFits : public DataIO
{
	private:

		pthread_mutex_t mutex;
		fitsfile *infile, *outfile;
		int fitsioStatus, anynull;
		string readKeywordStr(const char* keyword);
		long readKeywordLong(const char* keyword);
		float readKeywordFloat(const char* keyword);

		int nfields;
		float *x_phase_centre;
		float *y_phase_centre;

	public:
		DataIOFits(const char* infilename, const char* outfilename,
				pthread_mutex_t* mutex);
		~DataIOFits() {};
		int nvis();

		int readChunk(Chunk& chunk);
		void writeChunk(Chunk& chunk);

		int nPointings();
		float xPhaseCentre(int fieldID);
		float yPhaseCentre(int fieldID);
		void setPhaseCentre(int fieldID, double x, double y);
};

#endif // inclusion guard
