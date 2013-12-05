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
using casa::IPosition;

class Chunk;

class msio
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

	public:
		msio(const char* msinfile, const char* msoutfile, 
		          pthread_mutex_t* mutex);
		~msio();
		int nvis();

		int readChunk(Chunk& chunk);
		void writeChunk(Chunk& chunk);

		int nPointings();
		float xPhaseCentre(int id);
		float yPhaseCentre(int id);
		void setPhaseCenter(int fieldID, casa::Quantity x, casa::Quantity y);

};

#endif //inclusion guard
