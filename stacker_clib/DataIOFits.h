#include <pthread.h>
#include <fitsio.h>
#include <exception>
#include <stdexcept>

#include "DataIO.h"
using std::string;

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

		pthread_mutex_t* mutex;
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
		~DataIOFits();
		int nvis();

		int readChunk(Chunk& chunk);
		void writeChunk(Chunk& chunk);

		int nPointings();
		float xPhaseCentre(int fieldID);
		float yPhaseCentre(int fieldID);
		void setPhaseCentre(int fieldID, double x, double y);
};

#endif // inclusion guard
