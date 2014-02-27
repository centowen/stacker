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
		chunk.inVis[i].data.resize(1,1);
		chunk.inVis[i].weight.resize(1);
	}
	return chunk.size();
}

int msio::readChunkIteratorbased(Chunk& chunk)
{
	Vector<double> uvw;
	if(currentVisibility >= nvis()-1)
	{
		return 0;
	}
	else if(currentVisibility+chunk.size() > nvis())
		chunk.setSize(nvis()-currentVisibility);

	int uvrow;
	for(int i = 0; i < chunk.size(); i++)
	{
		uvrow = i+currentVisibility;
		chunk.inVis[i].index = uvrow;
	}
	for(int i = 0; i < chunk.size(); i++)
	{
		chunk.inVis[i].spw = msincols->dataDescId()(chunk.inVis[i].index);
	}
	for(int i = 0; i < chunk.size(); i++)
	{
		chunk.inVis[i].data = datainit->matrix();
		datainit->next();
	}

// 	casa::MatrixIterator<Complex>* datait = (casa::MatrixIterator<Complex>*) msc.data().getColumn().makeIterator(2);
// 		chunk.inVis[i].data = msincols.data()(uvrow);
	for(int i = 0; i < chunk.size(); i++)
	{
// 		uvw = msincols.uvw()(uvrow);
		uvw = uvwinit->vector();
		chunk.inVis[i].u = float(uvw[0]);
		chunk.inVis[i].v = float(uvw[1]);
		chunk.inVis[i].w = float(uvw[2]);
		uvwinit->next();
	}
	for(int i = 0; i < chunk.size(); i++)
	{
		chunk.inVis[i].weight  = weightinit->vector();
		weightinit->next();
	}
	for(int i = 0; i < chunk.size(); i++)
	{
		chunk.inVis[i].freq = freq[chunk.inVis[i].spw];
// 		chunk.inVis[i].freq  = new double[nchan];
// 		for(int col = 0; col < nchan; col++)
// 			chunk.inVis[i].freq[col] = freq[chunk.inVis[i].spw][col];
	}
	for(int i = 0; i < chunk.size(); i++)
	{
		uvrow = i+currentVisibility;
		chunk.inVis[i].fieldID = msincols->fieldId()(uvrow);
	}
	currentVisibility += chunk.size();
	return chunk.size();
}
int msio::readChunkSimple(Chunk& chunk)
{
	Vector<double> uvw;
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

	int uvrow = 0;
	for(int i = 0; i < chunk.size(); i++)
	{
		uvrow = i+currentVisibility;
		chunk.inVis[i].index = uvrow;
		chunk.inVis[i].data = msincols->data()(uvrow);
		uvw = msincols->uvw()(uvrow);
		chunk.inVis[i].u = float(uvw[0]);
		chunk.inVis[i].v = float(uvw[1]);
		chunk.inVis[i].w = float(uvw[2]);
		chunk.inVis[i].fieldID = msincols->fieldId()(uvrow);
		chunk.inVis[i].weight  = msincols->weight()(uvrow);
		chunk.inVis[i].data.unique();
		chunk.inVis[i].weight.unique();


		chunk.inVis[i].spw = msincols->dataDescId()(chunk.inVis[i].index);
		chunk.inVis[i].freq = freq[chunk.inVis[i].spw];
// 		chunk.inVis[i].freq  = new double[nchan];
// 		for(int col = 0; col < nchan; col++)
// 			chunk.inVis[i].freq[col] = freq[chunk.inVis[i].spw][col];
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
		msoutcols->data().put(chunk.outVis[i].index, chunk.outVis[i].data);
	}
	for(int i = 0; i < chunk.size(); i++)
	{
		msoutcols->weight().put(chunk.outVis[i].index, chunk.outVis[i].weight);
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

void msio::setPhaseCenter(int fieldID, casa::Quantity x, casa::Quantity y)
{
	casa::Array<double> newPhaseCentre(casa::IPosition(2,2,1));

	if(!msout)
		return;

	newPhaseCentre(IPosition(2,0,0)) = x.getValue("rad");
	newPhaseCentre(IPosition(2,1,0)) = y.getValue("rad");

	msoutcols->field().phaseDir().put(fieldID, newPhaseCentre);
	msoutcols->field().referenceDir().put(fieldID, newPhaseCentre);
	msoutcols->field().delayDir().put(fieldID, newPhaseCentre);
}

