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

#include <string>
#include <exception>
#include <stdexcept>

using std::string;
class Chunk;

#ifndef __DATA_H__
#define __DATA_H__

class fileException : public std::runtime_error {
	public:
		static const int OPEN = 0;
		static const int MOSAIC = 1;
		static const int HEADER_INFO_MISSING = 2;
		static const int CORRECTED_DATA_MISSING = 3;
		static const int MODEL_DATA_MISSING = 4;

	private:
		string message;
		int code;

	public:
		fileException(int code, string message) throw () : 
			std::runtime_error(message.c_str()), code(code){};
		~fileException() throw () {};

		int getCode() {return code; };
};

/* Base class for reading data from disk.
 *
 * DataIO is not thread-safe!
 */
class DataIO
{
	private:
		static int id_counter;
	public:
		DataIO();
		virtual ~DataIO();

		const int dataset_id;

		virtual size_t nvis() = 0;
		virtual int readChunk(Chunk& chunk) = 0;
		virtual void writeChunk(Chunk& chunk) = 0;
		virtual int nPointings() = 0;
		virtual float xPhaseCentre(int fieldID) = 0;
		virtual float yPhaseCentre(int fieldID) = 0;
		virtual void setPhaseCentre(int fieldID, double x, double y) = 0;

		virtual size_t nStokes() = 0;
		virtual size_t nChan() = 0;
		virtual size_t nSpw() = 0;
		virtual float* getFreq(int spw) = 0;

};

#endif //inclusion guard

