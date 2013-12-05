#include "MSComputer.h"
#include "Coords.h"
#include "PrimaryBeam.h"
#include "msio.h"

#ifndef __STACK_CHUNK_COMPUTER_H__
#define __STACK_CHUNK_COMPUTER_H__


class StackChunkComputer: public ChunkComputer
{
	private:
		Coords* coords;
		PrimaryBeam* pb;
		int stackingMode;
		bool redoWeights;

	public:
		void setStackingMode(int mode);
		StackChunkComputer(Coords* coords, PrimaryBeam* pb);

		// Called from computer and allows to access data,
		// unlike normal constructor which is called before computer
		// is created.
		void preCompute(msio* ms);
		virtual void computeChunk(Chunk* chunk);
		void postCompute(msio* ms);
};

#endif // inclusion guard

