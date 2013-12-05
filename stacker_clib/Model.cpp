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

using std::endl;
using std::cout;

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

void Model::compute(msio* ms, PrimaryBeam* pb)
{
	nPointings = (int)ms->nPointings();
	ComponentList cl(Path(clfile.c_str()));

    vector<double>* cx = new vector<double>[nPointings];
    vector<double>* cy = new vector<double>[nPointings];
    vector<double>* cflux = new vector<double>[nPointings];
    vector<double>* csize = new vector<double>[nPointings];

	nStackPoints = new int[nPointings];
	for(int i = 0; i < nPointings; i++)
	{
		nStackPoints[i] = 0;
	}

	double totFlux = 0.;
	std::cout << "Number of components in cl: " << cl.nelements() << std::endl;
	for(int i = 0; i < cl.nelements(); i++)
	{
		SkyComponent sc(cl.component(i));
		MDirection dir(sc.shape().refDirection());
		double x, y, flux, size;

		x = dir.getAngle().getValue("rad")[0];
		y = dir.getAngle().getValue("rad")[1];
		flux = double(sc.flux().value(casa::Stokes::I, false).getValue(Unit("Jy")));
		size = 0.;
		if(sc.shape().type() == 1)
			size = sc.shape().parameters()[0];
		totFlux += flux;

		for(int fieldID = 0; fieldID < nPointings; fieldID++)
		{
			double dx = 1.*((x - ms->xPhaseCentre(fieldID))*cos(y));
			double dy = 1.*(asin(sin(y)/cos(dx)) - ms->yPhaseCentre(fieldID));


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


	std::cout << "total flux: " << totFlux << std::endl;

	x = new double*[nPointings];
	y = new double*[nPointings];
	omega_x = new double*[nPointings];
	omega_y = new double*[nPointings];
	omega_z = new double*[nPointings];
	omega_size = new double*[nPointings];
	dx = new double*[nPointings];
	dy = new double*[nPointings];
	flux = new double*[nPointings];
	size = new double*[nPointings];

	for(int fieldID = 0; fieldID < nPointings; fieldID++)
		cout << "Number of model components in field " << fieldID << ": " << nStackPoints[fieldID] << endl;

	for(int fieldID = 0; fieldID < nPointings; fieldID++)
	{
		x[fieldID] = new double[nStackPoints[fieldID]];
		y[fieldID] = new double[nStackPoints[fieldID]];
		dx[fieldID] = new double[nStackPoints[fieldID]];
		dy[fieldID] = new double[nStackPoints[fieldID]];
		omega_x[fieldID] = new double[nStackPoints[fieldID]];
		omega_y[fieldID] = new double[nStackPoints[fieldID]];
		omega_z[fieldID] = new double[nStackPoints[fieldID]];
		omega_size[fieldID] = new double[nStackPoints[fieldID]];
		flux[fieldID] = new double[nStackPoints[fieldID]];
		size[fieldID] = new double[nStackPoints[fieldID]];

		for(int i = 0; i < nStackPoints[fieldID]; i++)
		{
			x[fieldID][i] = cx[fieldID][i];
			y[fieldID][i] = cy[fieldID][i];
			flux[fieldID][i] = cflux[fieldID][i];
			size[fieldID][i] = csize[fieldID][i];

			dx[fieldID][i] = (x[fieldID][i] - ms->xPhaseCentre(fieldID))*cos(y[fieldID][i]);
			dy[fieldID][i] = asin(sin(y[fieldID][i])/cos(dx[fieldID][i])) - ms->yPhaseCentre(fieldID);
			//          cout << "(dx,dy) (" << dx[fieldID][i] << ", " << dy[fieldID][i] << ") " << 180/pi*3600*sqrt(dx[fieldID][i]*dx[fieldID][i]+dy
			//          dx[fieldID][i] = (x[fieldID][i] - x_phase_centre[fieldID])*cos(y[fieldID][i]);
			//          dy[fieldID][i] = asin(sin(y[fieldID][i])) - y_phase_centre[fieldID];

			omega_x[fieldID][i] = 2*pi*sin(dx[fieldID][i])/casa::C::c;
			omega_y[fieldID][i] = 2*pi*sin(dy[fieldID][i])/casa::C::c;
			omega_z[fieldID][i] = 2*pi*(cos(sqrt(dx[fieldID][i]*dx[fieldID][i]+dy[fieldID][i]*dy[fieldID][i]))-1)/casa::C::c;
			omega_size[fieldID][i] = pow(pi*size[fieldID][i]/casa::C::c, 2) / (4 * log(2));
		}
	}

    delete[] cx;
    delete[] cy;
    delete[] cflux;
    delete[] csize;

	cout << "Model: " << endl;
	for(int fieldID = 0; fieldID < nPointings; fieldID++)
	{
		cout << "Field " << fieldID << ": " << endl;
		for(int i = 0; i < nStackPoints[fieldID]; i++)
		{
			cout << "omegas: (" 
			     << omega_x[fieldID][i] << ", " 
			     << omega_y[fieldID][i] << ", "
			     << omega_z[fieldID][i] << ", "
			     << omega_size[fieldID][i] << ")"
				 << " flux: " << flux[fieldID][i] 
				 << endl;
		}
	}
}
