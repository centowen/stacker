#include "PrimaryBeam.h"
#include <iostream>
#include "definitions.h"

#include <coordinates/Coordinates/CoordinateSystem.h>
#include <images/Images/ImageOpener.h>
#include <images/Images/ImageInfo.h>
#include <images/Images/ImageInterface.h>
#include <lattices/Lattices/LatticeBase.h>
#include <lattices/Lattices/Lattice.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>
#include <casa/Utilities/COWPtr.h>
#include <iostream>

using casa::IPosition;
using casa::ImageOpener;
using casa::CoordinateSystem;
using casa::ImageInfo;
using casa::COWPtr;
using casa::Slicer;
using casa::Array;
using casa::Vector;
// using casa::float;
using casa::Vector;
using casa::Quantum;
using casa::Unit;
using casa::uInt;

PrimaryBeam::PrimaryBeam(const char fileName[]) 
{
	ImageInterface<float>* interface = (ImageInterface<float>*)ImageOpener::openImage(fileName);
	cs = interface->coordinates();
	shape = interface->shape();
	nx = shape(0);
	ny = shape(1);

	IPosition start(4,0);
	IPosition dx(4,1,0,0,0);
	IPosition dy(4,0,1,0,0);

	const IPosition inputSliceShape(4,nx,1,1,1);
	Array<float> *buff_array = new Array<float>(inputSliceShape.nonDegenerate());
	Array<float> rowdata;
	COWPtr<Array<float> > inputArrPtr(buff_array);

	data = new float*[ny];

	for(uInt row = 0; row < ny; row++)
	{
		start(1) = row;
		interface->getSlice(inputArrPtr, Slicer(start, inputSliceShape), true);
		rowdata = *inputArrPtr;

		data[row] = new float[nx];

		for(uInt col = 0; col < nx; col++)
		{
			IPosition curpos = start+col*dx;
			data[row][col] = rowdata(curpos);
		}
	}

	x0 = float(cs.referenceValue()(0));
	y0 = float(cs.referenceValue()(1));

	if(cs.hasSpectralAxis())
		freq0 = float(cs.referenceValue()(cs.spectralAxisNumber()));
	else
		freq0 = float(0.);

	std::cout << "spectral: " << freq0 << std::endl;
	std::cout << "n: " << nx << ", " << ny << std::endl;

	px_x0 = float(cs.referencePixel()(0));
	px_y0 = float(cs.referencePixel()(1));
	this->dx = float(cs.increment()(0));
	this->dy = float(cs.increment()(1));

	delete interface;
}

PrimaryBeam::~PrimaryBeam()
{
	for(uInt row = 0; row < ny; row++)
		delete[] data[row];

	delete[] data;
}


float PrimaryBeam::calc(float x, float y, float freq)
{
	float freqcomp = 1.;
	if(freq > tol && freq0 > tol)
		freqcomp = freq0/freq;

	int px_x, px_y;
	px_x = int(freqcomp*(x-x0)/dx+px_x0);
	px_y = int(freqcomp*(y-y0)/dy+px_y0);

	if(px_x < 0 || px_x > nx) return 0.;
	if(px_y < 0 || px_y > ny) return 0.;
	return data[px_y][px_x];
// 	return data(IPosition(4,int(freqcomp*(x-x0)/dx+px_x0),int(freqcomp*(y-y0)/dy+px_y0),0,0));
}
