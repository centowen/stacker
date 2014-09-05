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
// Library to stack and modsub ms data.

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
