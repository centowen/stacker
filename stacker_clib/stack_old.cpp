// A simple test program to read in uv coverage and calculate maximum and miniumum baseline.


// includes/*{{{*/
#include <ms/MeasurementSets/MSTable.h>
#include <ms/MeasurementSets/MSColumns.h>
#include <ms/MeasurementSets/MeasurementSet.h>
#include <tables/Tables/TableDesc.h>
#include <casa/Arrays/VectorIter.h>
#include <casa/Arrays/MatrixIter.h>
#include <images/Images/ImageOpener.h>
#include <images/Images/ImageInfo.h>
#include <images/Images/ImageInterface.h>
#include <lattices/Lattices/Lattice.h>
#include <lattices/Lattices/LatticeBase.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <tables/Tables/ExprNodeSet.h>
#include <tables/Tables/ExprNode.h>
#include <casa/OS/Directory.h>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <pthread.h>
#include <queue>
#include <utility>
#include <string>
#include <string.h>
#include <stdio.h>
#include "PrimaryBeam.h"
#include "VisReader.h"
#include "Chunk.h"
#include "Coords.h"
#include "time.h"
#include "definitions.h"/*}}}*/

// using statements/*{{{*/
using std::cout;
using std::cerr;
using std::endl;
using std::fabs;
using std::exp;
using std::fstream;
using std::ifstream;
using std::getline;
using std::stringstream;
using std::vector;
using std::cos;
using std::sin;
using std::asin;
using std::queue;
using std::pair;
using std::string;

using casa::C::pi;

using casa::TableDesc;
using casa::IPosition;
using casa::Complex;
using casa::Float;
using casa::ImageInterface;
using casa::ImageOpener;
using casa::CoordinateSystem;
using casa::Directory;/*}}}*/

pthread_mutex_t mutex;
queue<int> chunksToWrite;
queue<int> chunksToCompute;
queue<int> freeChunks;
std::complex<float> sumVis;
float normSumVis;
int chunksDone, totalChunks;
Chunk** chunks;
bool allDataRead;

int stackingMode;
double beam;
bool weightByNvis;

PrimaryBeam* pb;
Coords* coords;

void computeChunk(Chunk* chunk) ;
Complex psf(double vis_real, double vis_image, double u, double v, double w, double dx, double dy, double freq);
void* computerThread(void* data)/*{{{*/
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
			return NULL;
		}

		pthread_mutex_unlock(&mutex);

		// Time to get to work if we found a chunk.
		if(chunkid>=0)
		{
			computeChunk(chunks[chunkid]);

			pthread_mutex_lock(&mutex);
			chunksToWrite.push(chunkid);
			pthread_mutex_unlock(&mutex);
		}

		// We could be down here because we are waiting for data to be read.
		// To be sure we don't use to much cpu while waiting lets relax a little.
		// This delay will not matter a lot even if there is data waiting.
		usleep(2000);
	}
}/*}}}*/

string to_string(int x)
{
	return dynamic_cast< std::ostringstream & >( \
		        ( std::ostringstream() << std::dec << x ) ).str();
};

int mainAction(queue<pair<int,string> > printQueue, VisReader& ms);
int mainActionSingleThread(queue<pair<int,string> > printQueue, VisReader& ms);
int mainActionMultiThread(queue<pair<int,string> > printQueue, VisReader& ms);

void calculate_chunk_average(Chunk& chunk);

double cpp_stack(const char* msinfile, const char* msoutfile, const char* pbfile, double* x, double* y, double* weight, int nstack);

// Container for data, only one container for library.
// If library is called multiple times in parallell this could be an issue.
__attribute__((constructor)) void initiate() {
	chunks = new Chunk*[N_CHUNK];
	for( int i = 0; i < N_CHUNK; i++)
	{
		chunks[i] = new Chunk(CHUNK_SIZE);
	}
};

// This is extern c to ensure no complications when called from outside library.
// Takes several input arguments:
//  - msinfile: The input ms file.
//  - msoutfile: The output ms file, can be the same as input ms file.
//  - pbfile: A casa image of the primary beam, used to calculate primary beam correction.
//  - x: x coordinate of each stacking position.
//  - y: y coordinate of each stacking position.
//  - weight: weight of each stacking position.
//  - nstack: length of x, y and weight lists.
//  - inStackingMode: Allows multiply each visibility with a psf function.
//  - inBeam: Parameter for inStackingMode psf function.
//  Returns average of all visibilities. Only useful is the case of simple data.
extern "C"{
	double stack(char* msinfile, char* msoutfile, char* pbfile, double* x, 
			     double* y, double* weight, int nstack, int inStackingMode, 
				 double inBeam)
	{
		char* l_msinfile;
		char* l_msoutfile; 
		char* l_coordfile; 
		double flux;
		int i;

		stackingMode = inStackingMode;
		if(stackingMode)
			beam = inBeam;
		else
			beam = 0.;

		flux = cpp_stack(msinfile, msoutfile, pbfile, x, y, weight, nstack);

		return flux;
	};
};

// Initiates actual stacking, called from the c version.
double cpp_stack(const char* msinfile, const char* msoutfile, const char* pbfile,
		double* x, double* y, double* weight, int nstack)/*{{{*/
{
	VisReader msdealer(msinfile, msoutfile, &mutex);
	pb = new PrimaryBeam(pbfile);
	coords = new Coords(msdealer, *pb, x, y, weight, nstack);

	if(msdealer.nPointings() > 1)
	{
		weightByNvis = true;
		cout << "Setting weightByNvis." << endl;
	}
	else
		weightByNvis = false;

	totalChunks = int(msdealer.nvis()/CHUNK_SIZE)+1;


	// Generate a queue of messages related to progress.
	// Specifies how many chunks needed to print a certain progress.
	queue<pair<int,string> > printQueue;
	for(int i = 0; i < 10; i++)
	{
		printQueue.push(pair<int,string>(int(i/10.*totalChunks), to_string(i*10)));
		for(int j = 0; j < 4; j++)
			printQueue.push(pair<int,string>(int(i/10.*totalChunks+(1+j)*totalChunks/50.), string(".")));
	}
	printQueue.push(pair<int,string>(totalChunks-1, "100%\n"));

// 	cout << "totchunk " << totalChunks << endl;
// 	cout << "Chunk size " << CHUNK_SIZE*sizeof(Visibility)*2./1024./1024. << "MiB" << endl;
// 	cout << "Number of visibilities " << msdealer.nvis() << endl;
// 	cout << "Time to start with stacking" << endl;

	sumVis = Complex(0.,0.);
	normSumVis = 0.;

	mainAction(printQueue, msdealer);

	msdealer.setPhaseCenter(0, 0., 0.);

	delete pb;
	delete coords;

	Complex averageFlux =  float(1. / normSumVis) * sumVis;

	return averageFlux.real();
}/*}}}*/

void computeChunk(Chunk* chunk) /*{{{*/
{
	for(int uvrow = 0; uvrow < chunk->size(); uvrow++)
	{
		// Shorthands to make code more readable.
		Visibility& inVis = chunk->inVis[uvrow];
		Visibility& outVis = chunk->outVis[uvrow];
		float &u = inVis.u;
		float &v = inVis.v;
		float &w = inVis.w;

		// Field ID is for which pointing the visibility is in.
		int fieldID = inVis.fieldID;

		// Data is in a matrix where columns are different frequencies
		// and rows are different polarizations.
		casa::Matrix<Complex> data = inVis.data;
		outVis.data = casa::Matrix<Complex>(inVis.data.shape());

		// For weights we only have a vector, ie. only one weight for all frequencies
		// and the rows represents different polarizations.
		casa::Vector<float> visWeight = inVis.weight;
		outVis.weight = casa::Vector<Float>(inVis.weight.shape());

		// Looping over frequency.
		for(int j = 0; j < inVis.data.ncolumn(); j++)
		{
			Complex dd(0,0);
			float dd_real = 0., dd_imag = 0.;
			float weightNorm = 0.;
			float phase;
			float freq = float(inVis.freq[j]);

			for(int i_p = 0; i_p < coords->nStackPoints[inVis.fieldID]; i_p++)
			{
				phase= -freq*(u*coords->omega_x[inVis.fieldID][i_p]+
				              v*coords->omega_y[inVis.fieldID][i_p]+
				              w*coords->omega_z[inVis.fieldID][i_p]);

				float weightbuff = coords->weight[fieldID][i_p];
				float pbcor = float(pb->calc(coords->dx[fieldID][i_p], coords->dy[fieldID][i_p], freq));

				dd_real = weightbuff*pbcor*cos(phase);
				dd_imag = weightbuff*pbcor*sin(phase);

				if(stackingMode == 0)
					dd += Complex(dd_real, dd_imag);
				else
					dd += psf(dd_real, dd_imag, u, v, w, coords->dx[fieldID][i_p], coords->dy[fieldID][i_p], freq);

				weightNorm += pbcor*pbcor*weightbuff;
			}

			if(weightNorm != 0)
			{
				dd /= weightNorm;
			}
			else
			{
				dd  = Complex(0,0);
			}

			// Looping over polarization.
			// dd does not need to be updated since it does not depend on polarization.
			for(int i = 0; i < inVis.data.nrow(); i++)
			{
				outVis.data(i,j) = dd*inVis.data(i,j);

				if(weightByNvis)
					if(weightNorm < 1e30)
						outVis.weight(i) = float(weightNorm)*inVis.weight(i);
					else
						outVis.weight(i) = float(0.0)*inVis.weight(i);
				else
					outVis.weight(i) = inVis.weight(i);
			}
		}

		outVis.fieldID = 0;
		outVis.index = inVis.index;
	}
}/*}}}*/


int mainAction(queue<pair<int,string> > printQueue, VisReader& ms)
{
// 	return mainActionSingleThread(printQueue, ms);
	return mainActionMultiThread(printQueue, ms);
}

int mainActionMultiThread(queue<pair<int,string> > printQueue, VisReader& ms)
{
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


	pthread_t threads[N_THREAD];
	pthread_mutex_init(&mutex, NULL);

	for(int i = 0; i < N_THREAD; i++)
	{
		pthread_create(&threads[i], NULL, computerThread, NULL);
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
			if(ms.readChunk(*chunks[chunkid]))
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
			calculate_chunk_average(*chunks[chunkid]);
			ms.writeChunk(*chunks[chunkid]);

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

	return 0;
}

// Single core version of mainActionMultiThread, clearly show what is done with
// the Chunks, read chunk, compute chunk, write chunk. Multi thread version 
// allows computeChunk to be run in parallell.
int mainActionSingleThread(queue<pair<int,string> > printQueue, VisReader& ms)
{
	Chunk chunk(CHUNK_SIZE);
	while(ms.readChunk(chunk))
	{
		computeChunk(&chunk);
		calculate_chunk_average(chunk);
	}
}

// Computes a running average to return at the end.
// Quick estimate for stacked flux as long as the source in a point source.
void calculate_chunk_average(Chunk& chunk)
{
	Complex sumVisLocal(float(0.),float(0.));
	float normSumVisLocal = 0.;
	for( int i = 0; i < chunk.size(); i++)
	{
		Complex aVis(0.,0.);
		for( int row = 0; row < chunk.outVis[i].data.nrow(); row ++)
			for( int column = 0; column < chunk.outVis[i].data.ncolumn(); column ++)
			{
				aVis += chunk.outVis[i].data(row, column);
			}
		aVis = float(1./(chunk.outVis[i].data.nrow()*chunk.outVis[i].data.ncolumn()))*aVis;
		sumVisLocal += float(chunk.outVis[i].weight(0))*aVis;
		normSumVisLocal += float(chunk.outVis[i].weight(0));
	}
	pthread_mutex_lock(&mutex);
	sumVis += sumVisLocal ;
	normSumVis += normSumVisLocal ;
	pthread_mutex_unlock(&mutex);
}

// Should indicate a psf to be used for stacking. Needs to be specfied in the
// uv plane. The simplest is of course dd = dd => real point. But this can 
// also be used to take things such as point in middle minus a circle around.
Complex psf(double vis_real, double vis_imag, double u, double v, double w, double dx, double dy, double freq)
{
	// Phi = -2*pi*freq/c*(u*dx+v*dy)
	// or more correct
	// Phi = -2*pi*freq/c*(u*sin(dx)+v*sin(dy)+(cos(sqrt(dx**2+dy**2))-1)*w)
	double Phi_u, Phi_v, Phi, Phi_v_h, Phi_u_h;

	if(stackingMode == 1)
	{
		Phi_u = -2*pi*freq/casa::C::c*u*beam/2.;
		Phi_v = -2*pi*freq/casa::C::c*v*beam/2.;
		Phi_u_h = -2*pi*freq/casa::C::c*(u+v)/sqrt(2)*beam/2.;
		Phi_v_h = -2*pi*freq/casa::C::c*(u-v)/sqrt(2)*beam/2.;
		vis_real = (1-(0.25*cos(Phi_u)+0.25*cos(Phi_v)+0.25*cos(Phi_u_h)+0.25*cos(Phi_v_h)))*vis_real;
	}
	else if(stackingMode == 2)
	{
		Phi = -2*pi*freq/casa::C::c/2.*(u*beam*cos(dx) + sin(dx)*beam*w);
		vis_real = (cos(Phi))*vis_real;
		vis_imag = (sin(Phi))*vis_imag;
// 		vis_real = (1-cos(Phi))*vis_real;
// 		vis_imag = (1-sin(Phi))*vis_real;
	}
	else if(stackingMode == 3)
	{
		Phi = 2*pi*freq/casa::C::c*beam*sqrt(u*u+v*v);
		vis_real = j0(Phi)*vis_real;
	}

	return Complex(vis_real, vis_imag);
}
