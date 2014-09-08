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
#include "Model.h"
#include <casa/Arrays/Array.h>
#include <images/Images/ImageInfo.h>
#include <images/Images/ImageInterface.h>
#include <images/Images/ImageOpener.h>
#include <lattices/Lattices/LatticeBase.h>
#include <lattices/Lattices/Lattice.h>
#include <casa/Arrays/Array.h>
#include <casa/Utilities/COWPtr.h>
#include <casa/OS/Directory.h>
#include <components/ComponentModels/ComponentList.h>
#include <components/ComponentModels/SkyComponent.h>
#include <components/ComponentModels/ComponentShape.h>
#include <components/ComponentModels/Flux.h>

using casa::Array;
using casa::Vector;
using casa::CoordinateSystem;
using casa::IPosition;
using casa::uInt;
using casa::COWPtr;
using casa::Slicer;
using casa::ImageInfo;
using casa::LatticeBase;
using casa::Lattice;
using casa::ImageOpener;
using casa::ImageInterface;
using casa::Path;
using casa::ComponentList;
using casa::SkyComponent;
using casa::ComponentShape;
using casa::Unit;
using casa::MDirection;

Model::Model(string file)
{
	clfile = file;

	nPointings = 0;
	nStackPoints = NULL;
	omega_x = NULL;
	omega_y = NULL;
	omega_z = NULL;
	omega_size = NULL;
	dx = NULL;
	dy = NULL;
	x = NULL;
	y = NULL;
	flux = NULL;
	size = NULL;
}

Model::~Model()
{
	if(nPointings > 0)
	{
		for(int i = 0; i < nPointings; i++)
		{
			delete[] omega_x[i];
			delete[] omega_y[i];
			delete[] omega_z[i];
			delete[] omega_size[i];
			delete[] dx[i];
			delete[] dy[i];
			delete[] x[i];
			delete[] y[i];
			delete[] flux[i];
			delete[] size[i];
		}
	}

	delete[] nStackPoints;
	delete[] omega_x;
	delete[] omega_y;
	delete[] omega_z;
	delete[] omega_size;
	delete[] dx;
	delete[] dy;
	delete[] x;
	delete[] y;
	delete[] flux;
	delete[] size;
}

void Model::compute(DataIO* ms, PrimaryBeam* pb)
{
	nPointings = (int)ms->nPointings();
	ComponentList cl(Path(clfile.c_str()));

    vector<float>* cx = new vector<float>[nPointings];
    vector<float>* cy = new vector<float>[nPointings];
    vector<float>* cflux = new vector<float>[nPointings];
    vector<float>* csize = new vector<float>[nPointings];

	nStackPoints = new int[nPointings];
	for(int i = 0; i < nPointings; i++)
	{
		nStackPoints[i] = 0;
	}

	double totFlux = 0.;

	for(int i = 0; i < cl.nelements(); i++)
	{
		SkyComponent sc(cl.component(i));
		MDirection dir(sc.shape().refDirection());
		float x, y, flux, size;

		x = float(dir.getAngle().getValue("rad")[0]);
		y = float(dir.getAngle().getValue("rad")[1]);
		flux = float(sc.flux().value(casa::Stokes::I, false).getValue(Unit("Jy")));
		size = float(0.);
		if(sc.shape().type() == 1)
			size = float(sc.shape().parameters()[0]);
		totFlux += flux;

		for(int fieldID = 0; fieldID < nPointings; fieldID++)
		{
			float dx = 1.*((x - float(ms->xPhaseCentre(fieldID)))*cos(y));
			float dy = 1.*(asin(sin(y)/cos(dx)) - float(ms->yPhaseCentre(fieldID)));


			if(pb->calc(dx,  dy )> 0.01)
			{
				cx[fieldID].push_back(x);
				cy[fieldID].push_back(y);
				cflux[fieldID].push_back(flux);
				csize[fieldID].push_back(size);
				nStackPoints[fieldID] ++;
			}
		}
	}


	x = new float*[nPointings];
	y = new float*[nPointings];
	omega_x = new float*[nPointings];
	omega_y = new float*[nPointings];
	omega_z = new float*[nPointings];
	omega_size = new float*[nPointings];
	dx = new float*[nPointings];
	dy = new float*[nPointings];
	flux = new float*[nPointings];
	size = new float*[nPointings];

	for(int fieldID = 0; fieldID < nPointings; fieldID++)
	{
		x[fieldID] = new float[nStackPoints[fieldID]];
		y[fieldID] = new float[nStackPoints[fieldID]];
		dx[fieldID] = new float[nStackPoints[fieldID]];
		dy[fieldID] = new float[nStackPoints[fieldID]];
		omega_x[fieldID] = new float[nStackPoints[fieldID]];
		omega_y[fieldID] = new float[nStackPoints[fieldID]];
		omega_z[fieldID] = new float[nStackPoints[fieldID]];
		omega_size[fieldID] = new float[nStackPoints[fieldID]];
		flux[fieldID] = new float[nStackPoints[fieldID]];
		size[fieldID] = new float[nStackPoints[fieldID]];

		for(int i = 0; i < nStackPoints[fieldID]; i++)
		{
			x[fieldID][i] = cx[fieldID][i];
			y[fieldID][i] = cy[fieldID][i];
			flux[fieldID][i] = cflux[fieldID][i];
			size[fieldID][i] = csize[fieldID][i];

			dx[fieldID][i] = (x[fieldID][i] - float(ms->xPhaseCentre(fieldID)))*cos(y[fieldID][i]);
			dy[fieldID][i] = asin(sin(y[fieldID][i])/cos(dx[fieldID][i])) - float(ms->yPhaseCentre(fieldID));

			omega_x[fieldID][i] = 2*pi*sin(dx[fieldID][i])/casa::C::c;
			omega_y[fieldID][i] = 2*pi*sin(dy[fieldID][i])/casa::C::c;
// 			omega_z[fieldID][i] = 2*pi*(cos(sqrt(dx[fieldID][i]*dx[fieldID][i]+dy[fieldID][i]*dy[fieldID][i]))-1)/casa::C::c;
			omega_z[fieldID][i] = 2*pi*(sqrt(1-dx[fieldID][i]*dx[fieldID][i]-dy[fieldID][i]*dy[fieldID][i])-1)/casa::C::c;
			omega_size[fieldID][i] = pow(pi*size[fieldID][i]/casa::C::c, 2) / (4 * log(2));
		}
	}

    delete[] cx;
    delete[] cy;
    delete[] cflux;
    delete[] csize;

}
