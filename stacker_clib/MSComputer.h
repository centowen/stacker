/***
 * MSComputer
 *
 * Runs a computer for all visbilities in a CASA measurement set. Will use
 * pthreads to use maximum number of cores. Will also handle all file io 
 * internally.
 ***/

#include <queue>
#include <string>
#include <pthread.h>

#include "definitions.h"
#include "msio.h"

#ifndef __MS_COMPUTER_H__
#define __MS_COMPUTER_H__

using std::queue;
using std::pair;
using std::string;

class Chunk;


class ChunkComputer
{
	public:
		virtual void computeChunk(Chunk* chunk) = 0;

		// Called before and after actual computation.
		// Allows full access to all data,
		virtual void preCompute(msio* ms) = 0;
		virtual void postCompute(msio* ms) = 0;
};

class MSComputer
{
	private:
		ChunkComputer* cc;
		Chunk** chunks;

		msio* ms;

		bool allDataRead;
		pthread_mutex_t mutex;

		queue<int> chunksToWrite, chunksToCompute, freeChunks;
		queue<pair<int,string> > printQueue;
		int chunksDone, totalChunks;
		string to_string(int x)
		{
			return dynamic_cast< std::ostringstream & >( \
						( std::ostringstream() << std::dec << x ) ).str();
		};

	public:
		MSComputer(ChunkComputer* cc, const char* msinfile, const char* msoutfile);
		~MSComputer();

		void run();
		static void* startComputerThread(void* data);
		void computerThread();
		msio* getMS();

};

#endif // inclusion guard
