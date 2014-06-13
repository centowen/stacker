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
// includes/*{{{*/
#include "definitions.h"
#include "MSComputer.h"
#include "Chunk.h"
#include <iostream>
/*}}}*/

using std::cout;
using std::endl;

MSComputer::MSComputer(ChunkComputer* cc, const char* infile, const char* outfile)/*{{{*/
{
	this->cc = cc;

	pthread_mutex_init(&mutex, NULL);

	chunks = new Chunk*[N_CHUNK];
	for( int i =0; i < N_CHUNK; i++)
		chunks[i] = new Chunk(CHUNK_SIZE);

	if(strcmp(".ms", infile+strlen(infile)-3) == 0
			||strcmp(".ms/", infile+strlen(infile)-4) == 0)
	{
		cout << "is ms" << endl;
		data = (DataIO*)(new msio(infile, outfile, &mutex));
	}
	else
		data = (DataIO*)(new DataIOFits(infile, outfile, &mutex));
}/*}}}*/

MSComputer::~MSComputer()/*{{{*/
{
	for( int i =0; i < N_CHUNK; i++)
		delete chunks[i];
	delete[] chunks;

	delete data;
}/*}}}*/

float MSComputer::run()/*{{{*/
{
	sumvis = 0.;
	nvis = 0.;

	totalChunks = int(data->nvis()/CHUNK_SIZE)+1;
	cout << "Total number of chunks should be " << totalChunks << endl;
	//
	// Generate a queue of messages related to progress.
	// Specifies how many chunks needed to print a certain progress.
	for(int i = 0; i < 10; i++)
	{
		printQueue.push(pair<int,string>(int(i/10.*totalChunks), to_string(i*10)));
		for(int j = 0; j < 4; j++)
			printQueue.push(pair<int,string>(int(i/10.*totalChunks+(1+j)*totalChunks/50.), string(".")));
	}
	printQueue.push(pair<int,string>(totalChunks-1, "100%\n"));

	// This variable ensures that threads are not closed
	// if there is a temporary hickup in data read.
	// Threads can only be terminated if all data is read,
	// even if that means that they have to sit idle for a while.
	// Also accessed from computer thread to find out if it is time to finish.
	allDataRead = false;

	// This counter is used for printing progress.
	// Not important for actual calculations.
	chunksDone = 0;

	// No chunks can be in use at start,
	// cleaning up from possible earlier run.
	while(!freeChunks.empty())
	{
		freeChunks.pop();
	}
	for( int i = 0; i < N_CHUNK; i++)
	{
		freeChunks.push(i);
	}

	cout << "Initiate cc" << endl;
	cc->preCompute(data);

	pthread_t threads[N_THREAD];

	cout << "Starting threads" << endl;
	for(int i = 0; i < N_THREAD; i++)
	{
		pthread_create(&threads[i], NULL, startComputerThread, (void*)this);
	}

	// Actual main job.
	// Chunks move through three queues.
	// 1. Free chunks are picked up by main thread and put into chunksToCompute
	// after reading data from disk into them.
	// 2. chunks from chunksToCompute are picked up by computer threads
	// and recalculated into stacked visibility before being put in chunksToWrite
	// 3. chunks from chunksToWrite are picked up by main thread and written to disk
	// then the chunks are put back into freeChunks
	//
	// All disk read and write is done by main thread.

	// Used to avoid stopping threads more than once and
	// also ensures that main thread waits for all computer threads
	// finish before moving on.
	bool threadsStillRunning = true;

	while(threadsStillRunning || !chunksToWrite.empty())
	{
		// Prints progress bar.
		while(!printQueue.empty() && chunksDone >= printQueue.front().first)
		{
			cout << printQueue.front().second << std::flush;
			printQueue.pop();
		}

		// Disk read and write.
		// Lock mutex to ensure queues don't change in the middle.
		// Would create complicated, possibly unpredictable behaviour.
		pthread_mutex_lock(&mutex);
		if(!allDataRead && !freeChunks.empty())
		{
			int chunkid = freeChunks.front();
			freeChunks.pop();
			pthread_mutex_unlock(&mutex);
			// Unlock mutex while reading from disk, 
			// chunk removed from queues ensure no one else
			// can access it. 
			if(data->readChunk(*chunks[chunkid]))
			{
				pthread_mutex_lock(&mutex);
				chunksToCompute.push(chunkid);

			}
			else
			{
				pthread_mutex_lock(&mutex);
				allDataRead = true;
			}
		}
		else if(!chunksToWrite.empty())
		{
			int chunkid = chunksToWrite.front();
			chunksToWrite.pop();
			pthread_mutex_unlock(&mutex);

			// Unlock mutex while writing to disk, 
			// chunk removed from queues ensure no one else
			// can access it.
// 			calculate_chunk_average(*chunks[chunkid]);
			data->writeChunk(*chunks[chunkid]);

			pthread_mutex_lock(&mutex);
			freeChunks.push(chunkid);
			chunksDone++;
		}
		pthread_mutex_unlock(&mutex);

		// When all data is read we only need to wait for computer thread to finish
		if(allDataRead && threadsStillRunning)
		{
			for(int i = 0; i < N_THREAD; i++)
			{
				pthread_join(threads[i], NULL);
			}
			threadsStillRunning = false;
		}

		usleep(2000);
	}

	cout << "clean up cc" << endl;
	cc->postCompute(data);

	cout << "sums: " << sumvis << ", " << nvis << endl;

	return sumvis/nvis;
}/*}}}*/

void* MSComputer::startComputerThread(void* computer)
{
	((MSComputer*)computer)->computerThread();
	return NULL;
}

void MSComputer::computerThread()/*{{{*/
{
	while(1)
	{
		// Any chunks to work on? Or is it time to die?
		int chunkid = -1;
		pthread_mutex_lock(&mutex);
		if(!chunksToCompute.empty())
		{
			chunkid = chunksToCompute.front();
			chunksToCompute.pop();
		}
		else if(allDataRead)
		{
			pthread_mutex_unlock(&mutex);
			return;
		}

		pthread_mutex_unlock(&mutex);

		float sumvis_chunk = 0.;
		// Time to get to work if we found a chunk.
		if(chunkid>=0)
		{
			sumvis_chunk = cc->computeChunk(chunks[chunkid]);

			pthread_mutex_lock(&mutex);
			if(abs(sumvis_chunk) > 0)
			{
				sumvis += sumvis_chunk;
				nvis += chunks[chunkid]->size();
			}
			chunksToWrite.push(chunkid);
			pthread_mutex_unlock(&mutex);
		}

		// We could be down here because we are waiting for data to be read.
		// To be sure we don't use to much cpu while waiting lets relax a little.
		// This delay will not matter a lot even if there is data waiting.
		usleep(2000);
	}
}/*}}}*/

DataIO* MSComputer::getMS()
{
	return data;
}

