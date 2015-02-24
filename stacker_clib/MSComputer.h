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

/***
 * MSComputer
 *
 * Runs a computer for all visbilities in a CASA measurement set. Will use
 * pthreads to use maximum number of cores. Will also handle all file io.
 ***/

#include <queue>
#include <string>
#include <pthread.h>

#include "definitions.h"
#include "DataIO.h"
#include "msio.h"
// #include "DataIOFits.h"

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
		virtual void preCompute(DataIO* data) = 0;
		virtual void postCompute(DataIO* data) = 0;
};

class MSComputer
{
	private:
		int n_thread_;
		ChunkComputer* cc;
		Chunk** chunks;

		DataIO* data;

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
		MSComputer(ChunkComputer* cc, 
				   int infiletype, const char* infilename, int infileoptions,
				   int outfiletype, const char* outfilename, int outfileoptions,
				   int n_thread = N_THREAD);
		~MSComputer();

		float run();
		static void* startComputerThread(void* data);
		void computerThread();
		DataIO* getMS();

};

#endif // inclusion guard
