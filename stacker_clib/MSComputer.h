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
#include "DataIO.h"
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
		virtual float computeChunk(Chunk* chunk) = 0;

		// Called before and after actual computation.
		// Allows full access to all data,
		virtual void preCompute(DataIO* data) = 0;
		virtual void postCompute(DataIO* data) = 0;
};

class MSComputer
{
	private:
		ChunkComputer* cc;
		Chunk** chunks;

		DataIO* data;

		bool allDataRead;
		pthread_mutex_t mutex;

		queue<int> chunksToWrite, chunksToCompute, freeChunks;
		queue<pair<int,string> > printQueue;
		int chunksDone, totalChunks;

		float sumvis;
		float nvis;

		string to_string(int x)
		{
			return dynamic_cast< std::ostringstream & >( \
						( std::ostringstream() << std::dec << x ) ).str();
		};

	public:
		MSComputer(ChunkComputer* cc, const char* msinfile, const char* msoutfile);
		~MSComputer();

		float run();
		static void* startComputerThread(void* data);
		void computerThread();
		DataIO* getMS();

};

#endif // inclusion guard
