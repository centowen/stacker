#include "MSComputer.h"
#include "Model.h"
#include "PrimaryBeam.h"
#include "msio.h"

#ifndef __MODSUB_CHUNK_COMPUTER_H__
#define __MODSUB_CHUNK_COMPUTER_H__

class ModsubChunkComputer: public ChunkComputer
{
	private:
		Model* model;
		PrimaryBeam* pb;

	public:
		ModsubChunkComputer(Model* model, PrimaryBeam* pb);
		~ModsubChunkComputer();

		// Called from computer and allows to access data,
		// unlike normal constructor which is called before computer
		// is created.
		void preCompute(msio* ms);
		virtual float computeChunk(Chunk* chunk);
		void postCompute(msio* ms);
};

#endif // inclusion guard

