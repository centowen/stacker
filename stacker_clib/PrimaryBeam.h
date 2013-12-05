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
// 		ImageInterface<float>* interface;
// 		Array<float> data, *buff_array;
		float** data;
		CoordinateSystem cs;
		IPosition shape;
		static const double pi = 3.141592653589793238462;
		int nx, ny;
		float x0, y0, dx, dy, px_x0, px_y0;
		float freq0;

	public:
		PrimaryBeam(const char fileName[]);
		~PrimaryBeam();
// 		void load(const char fileName[]);
		float calc(float x, float y, float freq = 0.);
};

#endif
