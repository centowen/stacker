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
#include "msio.h"
#include "Chunk.h"
#include <iostream>
#include <stdlib.h>

msio::msio(const char* msinfile, const char * msoutfile,
		             pthread_mutex_t* mutex )
{
	msin = new MeasurementSet(msinfile);
	msincols = new ROMSColumns(*msin);

	if(strlen(msoutfile) > 0)
	{
		msout = new MeasurementSet(msoutfile, casa::Table::Update);
		msoutcols = new MSColumns(*msout);
	}
	else
	{
		msout = NULL;
		msoutcols = NULL;
	}
	this->mutex = *mutex;
	currentVisibility = 0;
// 	datainit = (casa::MatrixIterator<Complex>*) msincols->data().getColumn().makeIterator(2);
// 	uvwinit = (casa::VectorIterator<double>*) msincols->uvw().getColumn().makeIterator(1);
// 	weightinit = (casa::VectorIterator<float>*) msincols->weight().getColumn().makeIterator(1);

	VectorIterator<double>* freqinit = (casa::VectorIterator<double>*) msincols->spectralWindow().chanFreq().getColumn().makeIterator(1);
	
	nspw = msincols->spectralWindow().nrow();
	freq  = new double*[nspw];
	for(int row =0 ; row < nspw; row++)
	{
		casa::Vector<double> freqbuff = freqinit->vector();
		nchan = freqbuff.shape()(0);
		freq[row] = new double[nchan];
		for(int col = 0; col < nchan; col++)
			freq[row][col] = freqbuff(col);
		freqinit->next();
	}

	nfields = msin->field().nrow();
	x_phase_centre = new float[nfields];
	y_phase_centre = new float[nfields];
	for(int fieldID = 0; fieldID < nfields; fieldID++)
	{
		Array<double> phase_centre = msincols->field().phaseDir()(fieldID);
		x_phase_centre[fieldID] = float(phase_centre(IPosition(2,0,0)));
		y_phase_centre[fieldID] = float(phase_centre(IPosition(2,1,0)));
	}

}

msio::~msio()
{
// 	msin.flush();
	delete msincols;
	delete msin;
	if(msout)
		delete msout, msoutcols;
// 	for(int i = 0; i < nspw; i++)
// 		delete freq[i];
// 	delete freq;
}

int msio::nvis()
{
	return msincols->data().nrow();
}

int msio::readChunk(Chunk& chunk)
{
// 	readChunkIteratorbased(chunk);
	readChunkSimple(chunk);
}

int msio::readChunkDummy(Chunk& chunk)
{
	if(currentVisibility >= nvis()-1)
		return 0;
	currentVisibility += chunk.size();
	int uvrow = 0;
	for(int i = 0; i < chunk.size(); i++)
	{
		uvrow = i+currentVisibility;
		chunk.inVis[i].index = uvrow;
		chunk.inVis[i].u = float(0.);
		chunk.inVis[i].v = float(0.);
		chunk.inVis[i].w = float(0.);
		chunk.inVis[i].fieldID = 0;
		chunk.inVis[i].spw = 0;
		chunk.inVis[i].freq = freq[chunk.inVis[i].spw];
	}
	return chunk.size();
}

int msio::readChunkIteratorbased(Chunk& chunk)
{
	readChunkDummy(chunk);
}
// {/*{{{*/
// 	Vector<double> uvw;
// 	if(currentVisibility >= nvis()-1)
// 	{
// 		return 0;
// 	}
// 	else if(currentVisibility+chunk.size() > nvis())
// 		chunk.setSize(nvis()-currentVisibility);
// 
// 	int uvrow;
// 	for(int i = 0; i < chunk.size(); i++)
// 	{
// 		uvrow = i+currentVisibility;
// 		chunk.inVis[i].index = uvrow;
// 	}
// 	for(int i = 0; i < chunk.size(); i++)
// 	{
// 		chunk.inVis[i].spw = msincols->dataDescId()(chunk.inVis[i].index);
// 	}
// 	for(int i = 0; i < chunk.size(); i++)
// 	{
// 		chunk.inVis[i].data = datainit->matrix();
// 		datainit->next();
// 	}
// 
// // 	casa::MatrixIterator<Complex>* datait = (casa::MatrixIterator<Complex>*) msc.data().getColumn().makeIterator(2);
// // 		chunk.inVis[i].data = msincols.data()(uvrow);
// 	for(int i = 0; i < chunk.size(); i++)
// 	{
// // 		uvw = msincols.uvw()(uvrow);
// 		uvw = uvwinit->vector();
// 		chunk.inVis[i].u = float(uvw[0]);
// 		chunk.inVis[i].v = float(uvw[1]);
// 		chunk.inVis[i].w = float(uvw[2]);
// 		uvwinit->next();
// 	}
// 	for(int i = 0; i < chunk.size(); i++)
// 	{
// 		chunk.inVis[i].weight  = weightinit->vector();
// 		weightinit->next();
// 	}
// 	for(int i = 0; i < chunk.size(); i++)
// 	{
// 		chunk.inVis[i].freq = freq[chunk.inVis[i].spw];
// // 		chunk.inVis[i].freq  = new double[nchan];
// // 		for(int col = 0; col < nchan; col++)
// // 			chunk.inVis[i].freq[col] = freq[chunk.inVis[i].spw][col];
// 	}
// 	for(int i = 0; i < chunk.size(); i++)
// 	{
// 		uvrow = i+currentVisibility;
// 		chunk.inVis[i].fieldID = msincols->fieldId()(uvrow);
// 	}
// 	currentVisibility += chunk.size();
// 	return chunk.size();
// }/*}}}*/
//
int msio::readChunkSimple(Chunk& chunk)
{
	Vector<double> uvw;
	Matrix<Complex> data;
	Vector<float> weight;
	pthread_mutex_lock(&mutex);
	int currentVisibility = this->currentVisibility;
	if(currentVisibility >= nvis())
		return 0;
	else if(currentVisibility+chunk.size() > nvis())
	{
		chunk.setSize(nvis()-currentVisibility);
	}
	this->currentVisibility += chunk.size();
	pthread_mutex_unlock(&mutex);

	int uvrow = 0, nchan, nstokes;
	for(int i = 0; i < chunk.size(); i++)
	{
		uvrow = i+currentVisibility;
		data = msincols->data()(uvrow);
		weight = msincols->weight()(uvrow);
		uvw = msincols->uvw()(uvrow);

		chunk.inVis[i].index = uvrow;

        nchan = data.ncolumn();
        nstokes = data.nrow();
        if(nchan != chunk.inVis[i].nchan || nstokes != chunk.inVis[i].nstokes)
        {
            chunk.inVis[i].nchan = nchan;
            chunk.inVis[i].nstokes = nstokes;
            chunk.inVis[i].data_real = new float[nchan*nstokes];
            chunk.inVis[i].data_imag = new float[nchan*nstokes];
            chunk.inVis[i].weight = new float[nstokes];
        }
        for(int stokes = 0; stokes < nstokes; stokes++)
        {
            chunk.inVis[i].weight[stokes] = float(weight(stokes));
            for(int chan = 0; chan < nchan; chan++)
            {
                chunk.inVis[i].data_real[nchan*stokes+chan] = float(std::real(data(stokes,chan)));
                chunk.inVis[i].data_imag[nchan*stokes+chan] = float(std::imag(data(stokes,chan)));
            }
        }

		chunk.inVis[i].u = float(uvw[0]);
		chunk.inVis[i].v = float(uvw[1]);
		chunk.inVis[i].w = float(uvw[2]);
		chunk.inVis[i].fieldID = msincols->fieldId()(uvrow);

		chunk.inVis[i].spw = msincols->dataDescId()(chunk.inVis[i].index);
		chunk.inVis[i].freq = freq[chunk.inVis[i].spw];

        if(nchan != chunk.outVis[i].nchan || nstokes != chunk.outVis[i].nstokes)
        {
            chunk.outVis[i].nchan = nchan;
            chunk.outVis[i].nstokes = nstokes;
            chunk.outVis[i].data_real = new float[nchan*nstokes];
            chunk.outVis[i].data_imag = new float[nchan*nstokes];
            chunk.outVis[i].weight = new float[nstokes];
        }
	}
	return chunk.size();
}

void msio::writeChunk(Chunk& chunk)
{
	if(msout == NULL)
		return;

	Vector<double> uvw;
	for(int i = 0; i < chunk.size(); i++)
	{
		int nchan = chunk.outVis[i].nchan, 
			nstokes = chunk.outVis[i].nstokes;
		Matrix<Complex> data(nstokes, nchan);
		for(int chan = 0; chan < nchan; chan++)
		{
			for(int stokes = 0; stokes < nstokes; stokes++)
			{
				Complex vis = Complex(chunk.outVis[i].data_real[stokes*nchan+chan],
						              chunk.outVis[i].data_imag[stokes*nchan+chan]);
				data(stokes, chan) = vis;
			}
		}
		msoutcols->data().put(chunk.outVis[i].index, data);
	}
	for(int i = 0; i < chunk.size(); i++)
	{
		int nstokes = chunk.outVis[i].nstokes;
		Vector<Float> weight(nstokes);
		for(int stokes = 0; stokes<nstokes; stokes++)
			weight(stokes) = chunk.outVis[i].weight[stokes];

		msoutcols->weight().put(chunk.outVis[i].index, weight);
	}
	for(int i = 0; i < chunk.size(); i++)
	{
		msoutcols->fieldId().put(chunk.outVis[i].index, chunk.outVis[i].fieldID);
	}
}

int msio::nPointings()
{
	return nfields;
}

float msio::xPhaseCentre(int id)
{
	return x_phase_centre[id];
}

float msio::yPhaseCentre(int id)
{
	return y_phase_centre[id];
}

void msio::setPhaseCentre(int fieldID, double x, double y)
{
	casa::Array<double> newPhaseCentre(casa::IPosition(2,2,1));

	if(!msout)
		return;

	newPhaseCentre(IPosition(2,0,0)) = x;
	newPhaseCentre(IPosition(2,1,0)) = y;

	msoutcols->field().phaseDir().put(fieldID, newPhaseCentre);
	msoutcols->field().referenceDir().put(fieldID, newPhaseCentre);
	msoutcols->field().delayDir().put(fieldID, newPhaseCentre);
}

