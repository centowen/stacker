#ifndef __PRIMARYBEAM_H__
#define __PRIMARYBEAM_H__

class PrimaryBeam
{
	private:
		static const double pi = 3.141592653589793238462;

	public:
		PrimaryBeam();
		~PrimaryBeam();
		virtual float calc(float x, float y, float freq = 0.) = 0;
};

class ConstantPrimaryBeam: public PrimaryBeam
{
	public:
		ConstantPrimaryBeam();
		~ConstantPrimaryBeam();
		virtual float calc(float x, float y, float freq = 0.);
};

class ImagePrimaryBeam: public PrimaryBeam
{
	protected:
		float** data;
		int nx, ny;
		float x0, y0, dx, dy, px_x0, px_y0;
		float freq0;

	public:
		ImagePrimaryBeam(const char fileName[]);
		~ImagePrimaryBeam();

		virtual float calc(float x, float y, float freq = 0.);
};

#endif
