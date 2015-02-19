#include <vector>

#include "DataIO.h"
#include "Chunk.h"

#ifndef __CACHED_DATAIO_H__
#define __CACHED_DATAIO_H__

using std::vector;

class CachedDataIO : public DataIO
{
	private:
		size_t max_chunks_;
		vector<Chunk> cache;
		vector<Chunk>::iterator cache_iterator;

		size_t nchan, nspw, nstokes;
		int n_vis;
		int nfields;
		float* x_phase_centre;
		float* y_phase_centre;
		float* freq;

	public:
		CachedDataIO(DataIO* dataio, size_t max_chunks);
		int nvis();

		int readChunk(Chunk& chunk);
		void writeChunk(Chunk& chunk);

		int nPointings();
		float xPhaseCentre(int fieldID);
		float yPhaseCentre(int fieldID);
		void setPhaseCentre(int fieldID, double x, double y);

		size_t nChan();
		size_t nSpw();
		size_t nStokes();
		float* getFreq(int spw);
		void restart();
};

#endif // inclusion guard
