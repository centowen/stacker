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
using casa::Double;
using casa::Vector;
using casa::Quantum;
using casa::Unit;

PrimaryBeam::PrimaryBeam(const char fileName[]) 
{
	interface = (ImageInterface<float>*)ImageOpener::openImage(fileName);
	cs = interface->coordinates();
	shape = interface->shape();
	const IPosition inputSliceShape(4,shape(0),shape(1),1,1);
	buff_array = new Array<float>(inputSliceShape.nonDegenerate());
	COWPtr<Array<float> > inputArrPtr(buff_array);
	IPosition start(4,0);

	interface->getSlice(inputArrPtr, Slicer(start, inputSliceShape), true);
	data = *inputArrPtr;

	x0 = cs.referenceValue()(0);
	y0 = cs.referenceValue()(1);

	if(cs.hasSpectralAxis())
		freq0 = cs.referenceValue()(cs.spectralAxisNumber());
	else
		freq0 = 0.;

	std::cout << "spectral: " << freq0 << std::endl;

	px_x0 = cs.referencePixel()(0);
	px_y0 = cs.referencePixel()(1);
	dx = cs.increment()(0);
	dy = cs.increment()(1);

}

PrimaryBeam::~PrimaryBeam()
{
	delete interface;
}

double PrimaryBeam::calc(double x, double y, double freq)
{
	double freqcomp = 1.;
	if(freq > tol && freq0 > tol)
		freqcomp = freq0/freq;

	if(int(freqcomp*(x-x0)/dx+px_x0) < 0 || int(freqcomp*(x-x0)/dx+px_x0) > shape(0)) return 0.;
	if(int(freqcomp*(y-y0)/dy+px_y0) < 0 || int(freqcomp*(y-y0)/dy+px_y0) > shape(1)) return 0.;
	return data(IPosition(4,int(freqcomp*(x-x0)/dx+px_x0),int(freqcomp*(y-y0)/dy+px_y0),0,0));
}
