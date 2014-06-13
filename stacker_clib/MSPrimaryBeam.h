#include "PrimaryBeam.h"

#include <casa/Arrays/Array.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <images/Images/ImageInterface.h>

#ifndef __MSPRIMARYBEAM_H__
#define __MSPRIMARYBEAM_H__

using casa::Array;
using casa::CoordinateSystem;
using casa::IPosition;
using casa::ImageInterface;

class MSPrimaryBeam: public ImagePrimaryBeam
{
	private:
		CoordinateSystem cs;

	public:
		MSPrimaryBeam(const char fileName[]);
		~MSPrimaryBeam();
};

#endif // End of inclusion guard
