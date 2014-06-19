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
#include <iostream>

#include "StackChunkComputer.h"
#include "Chunk.h"
#include "Coords.h"
#include "PrimaryBeam.h"

using std::cout;
using std::real;

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

float StackChunkComputer::computeChunk(Chunk* chunk) /*{{{*/
{
	float sum = 0., normsum = 0.;
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
		float* data_real = inVis.data_real;
		float* data_imag = inVis.data_imag;

		// For weights we only have a vector, ie. only one weight for all frequencies
		// and the rows represents different polarizations.
		float* visWeight = inVis.weight;
// 		outVis.weight = casa::Vector<Float>(inVis.weight.shape());

		// Looping over frequency.
		for(int j = 0; j < inVis.nchan; j++)
		{
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

				dd_real += weightbuff*pbcor*cos(phase);
				dd_imag += weightbuff*pbcor*sin(phase);

// 				dd += Complex(dd_real, dd_imag);
// 				if(stackingMode == 0)
// 					dd += Complex(dd_real, dd_imag);
// 				else
// 					dd += psf(dd_real, dd_imag, u, v, w, coords->dx[fieldID][i_p], coords->dy[fieldID][i_p], freq);

				weightNorm += pbcor*pbcor*weightbuff;
			}

			if(weightNorm != 0)
			{
				dd_real /= weightNorm;
				dd_imag /= weightNorm;
			}
			else
			{
				dd_real  = 0.;
				dd_imag  = 0.;
			}

			// Looping over polarization.
			// dd does not need to be updated since it does not depend on polarization.
			for(int i = 0; i < inVis.nstokes; i++)
			{
// 				outVis.data_real[i*outVis.nchan+j] = inVis.data_real[i*inVis.nchan+j];
// 				outVis.data_imag[i*outVis.nchan+j] = inVis.data_imag[i*inVis.nchan+j];
				outVis.data_real[i*outVis.nchan+j] = dd_real*inVis.data_real[i*inVis.nchan+j]
					                               - dd_imag*inVis.data_imag[i*inVis.nchan+j];
				outVis.data_imag[i*outVis.nchan+j] = dd_real*inVis.data_imag[i*inVis.nchan+j]
					                               + dd_imag*inVis.data_real[i*inVis.nchan+j];


				if(redoWeights)
					if(weightNorm < 1e30)
						outVis.weight[i] = float(weightNorm)*inVis.weight[i];
					else
						outVis.weight[i] = float(0.0)*inVis.weight[i];
				else
					outVis.weight[i] = inVis.weight[i];

				sum += outVis.data_real[i*inVis.nchan+j]*outVis.weight[i];
				normsum += outVis.weight[i];
			}
		}

		if(coords->nStackPoints[inVis.fieldID] > 0)
			outVis.fieldID = 0;
		else
			outVis.fieldID = 1;
		outVis.index = inVis.index;
	}
	float retval = 0;
	if( normsum > 0)
		retval = sum/normsum*chunk->size();
	return retval;
// 	return 1.;
}/*}}}*/

void StackChunkComputer::preCompute(DataIO* data)
{
	// If we only have one field in the data the weights do not need to be 
	// updated. 
	if(data->nPointings() > 1)
	{
		redoWeights = true;
	}

	coords->computeCoords(data, *pb);
}

void StackChunkComputer::postCompute(DataIO* data)
{
	// After stacking we have all visibilities in field 0.
	// Set centre to 0., 0. to indicate that coordinates at this point are
	// arbitrary.
	data->setPhaseCentre(0, 0., 0.);

	// Maybe it would be nice to also remove all the other fields?
	// Also should properly flag visibilities that got bad in stacking.
}
