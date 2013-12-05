#include <casa/Arrays/Matrix.h>
#include <iostream>

#include "StackChunkComputer.h"
#include "Chunk.h"
#include "Coords.h"
#include "PrimaryBeam.h"

using casa::Matrix;
using casa::Complex;
using casa::Float;
using std::cout;

StackChunkComputer::StackChunkComputer(Coords* coords, PrimaryBeam* pb)
{
	this->coords = coords;
	this->pb = pb;
	stackingMode = 0;
}

void StackChunkComputer::setStackingMode(int mode)
{
	stackingMode = mode;
}

void StackChunkComputer::computeChunk(Chunk* chunk) /*{{{*/
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
		Matrix<Complex> data = inVis.data;
		outVis.data = Matrix<Complex>(inVis.data.shape());

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

				dd += Complex(dd_real, dd_imag);
// 				if(stackingMode == 0)
// 					dd += Complex(dd_real, dd_imag);
// 				else
// 					dd += psf(dd_real, dd_imag, u, v, w, coords->dx[fieldID][i_p], coords->dy[fieldID][i_p], freq);

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

				if(redoWeights)
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

void StackChunkComputer::preCompute(msio* ms)
{
	// If we only have one field in the data the weights do not need to be 
	// updated. 
	if(ms->nPointings() > 1)
	{
		redoWeights = true;
	}

	coords->computeCoords(ms, *pb);
}

void StackChunkComputer::postCompute(msio* ms)
{
	// After stacking we have all visibilities in field 0.
	// Set centre to 0., 0. to indicate that coordinates at this point are
	// arbitrary.
	ms->setPhaseCenter(0, 0., 0.);

	// Maybe it would be nice to also remove all the other fields?
	// Also should properly flag visibilities that got bad in stacking.
}
