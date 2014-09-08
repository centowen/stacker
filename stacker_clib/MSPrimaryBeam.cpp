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
#include "MSPrimaryBeam.h"

#include <coordinates/Coordinates/CoordinateSystem.h>
#include <images/Images/ImageOpener.h>
#include <images/Images/ImageInfo.h>
#include <images/Images/ImageInterface.h>
#include <lattices/Lattices/LatticeBase.h>
#include <lattices/Lattices/Lattice.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>
#include <casa/Utilities/COWPtr.h>

using casa::IPosition;
using casa::ImageOpener;
using casa::CoordinateSystem;
using casa::ImageInfo;
using casa::COWPtr;
using casa::Slicer;
using casa::Array;
using casa::Vector;
using casa::Vector;
using casa::Quantum;
using casa::Unit;
using casa::uInt;

MSPrimaryBeam::MSPrimaryBeam(const char fileName[]) /*{{{*/
	: ImagePrimaryBeam(fileName)
{
	ImageInterface<float>* interface = (ImageInterface<float>*)ImageOpener::openImage(fileName);
	cs = interface->coordinates();
	IPosition shape = interface->shape();
	this->nx = shape(0);
	this->ny = shape(1);

	IPosition start(4,0);
	IPosition dx(4,1,0,0,0);
	IPosition dy(4,0,1,0,0);

	const IPosition inputSliceShape(4,nx,1,1,1);
	Array<float> *buff_array = new Array<float>(inputSliceShape.nonDegenerate());
	Array<float> rowdata;
	COWPtr<Array<float> > inputArrPtr(buff_array);

	this->data = new float*[ny];

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

	px_x0 = float(cs.referencePixel()(0));
	px_y0 = float(cs.referencePixel()(1));
	this->dx = float(cs.increment()(0));
	this->dy = float(cs.increment()(1));

	delete interface;
}/*}}}*/

MSPrimaryBeam::~MSPrimaryBeam()
{
	for(uInt row = 0; row < ny; row++)
		delete[] data[row];

	delete[] data;
}
