#include <casa/Arrays/Array.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <images/Images/ImageInterface.h>

#ifndef __PRIMARYBEAM_H__
#define __PRIMARYBEAM_H__

using casa::Array;
using casa::CoordinateSystem;
using casa::IPosition;
using casa::ImageInterface;

class PrimaryBeam
{
	private:
		static const double pi = 3.141592653589793238462;

	public:
		virtual float calc(float x, float y, float freq = 0.) = 0;
};

class ImagePrimaryBeam: public PrimaryBeam
{
	private:
		float** data;
		CoordinateSystem cs;
		IPosition shape;
		int nx, ny;
		float x0, y0, dx, dy, px_x0, px_y0;
		float freq0;

	public:
		ImagePrimaryBeam(const char fileName[]);
		~ImagePrimaryBeam();

		float calc(float x, float y, float freq = 0.);
};

#endif
