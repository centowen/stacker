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
#include <casa/Arrays/Matrix.h>
#include <iostream>

#include "ModsubChunkComputer.h"
#include "Chunk.h"
#include "Coords.h"
#include "PrimaryBeam.h"

using casa::Matrix;
using casa::Vector;
using casa::Complex;
using casa::Float;
using std::cout;


ModsubChunkComputer::ModsubChunkComputer(Model* model, PrimaryBeam* pb)
{
	this->model = model;
	this->pb = pb;
}

ModsubChunkComputer::~ModsubChunkComputer() {}

float ModsubChunkComputer::computeChunk(Chunk* chunk) /*{{{*/
{
	for(int uvrow = 0; uvrow < chunk->size(); uvrow++)
	{
// 		int fieldID = msc.fieldId()(uvrow);

		Visibility& inVis = chunk->inVis[uvrow];
		Visibility& outVis = chunk->outVis[uvrow];

		Matrix<Complex> data = inVis.data;
		Vector<float> visWeight = inVis.weight;
		int fieldID = inVis.fieldID;

		outVis.data = Matrix<Complex>(inVis.data.shape());
		outVis.weight = Vector<Float>(inVis.weight.shape());


		float &u = inVis.u;
		float &v = inVis.v;
		float &w = inVis.w;

		for(int j = 0; j < inVis.data.ncolumn(); j++)
		{
			for(int i = 0; i < inVis.data.nrow(); i++)
			{
				Complex dd(0,0);
				float dd_real = 0., dd_imag = 0.;
				float weightNorm = 0.;
				float d;
				float freq = float(inVis.freq[j]);

				for(int i_p = 0; i_p < model->nStackPoints[inVis.fieldID]; i_p++)
				{
					float extent = 1.;

// 					d = - freq*(u*omega_x[fieldID][i_p]+v*omega_y[fieldID][i_p]);
					d = freq*(u*model->omega_x[inVis.fieldID][i_p]+
						v*model->omega_y[inVis.fieldID][i_p]+
						w*model->omega_z[inVis.fieldID][i_p]);
					if( model->size[inVis.fieldID][i_p] > 1e-10)
						extent = exp(-freq*freq*(u*u + v*v)*model->omega_size[inVis.fieldID][i_p]);

					float pbcor = float(pb->calc(model->dx[fieldID][i_p], 
								                 model->dy[fieldID][i_p], 
												 freq));
 					dd_real += pbcor*model->flux[fieldID][i_p]*extent*cos(d);
 					dd_imag += pbcor*model->flux[fieldID][i_p]*extent*sin(d);

				}

				dd = Complex(dd_real, dd_imag);
				outVis.data(i,j) = inVis.data(i,j)-dd;
			}
		}

		outVis.weight = inVis.weight;
		outVis.data.unique();
		outVis.weight.unique();
		outVis.fieldID = inVis.fieldID;
		outVis.index = inVis.index;

//  		msc_out.data().put(uvrow, (casa::Array<Complex>)data);
//  		msc_out.weight().put(uvrow, (casa::Array<Float>)visWeight);
//  		msc_out.fieldId().put(uvrow, 0);
		
// 		msc_out.correctedData().put(i, (casa::Array<Complex>)data);

// 		it->next();
// 		visWeightit->next();
// 		datait->next();
	}
	return 0.;
}/*}}}*/

void ModsubChunkComputer::preCompute(DataIO* ms)
{
	model->compute(ms, pb);
}

void ModsubChunkComputer::postCompute(DataIO* ms)
{
}
