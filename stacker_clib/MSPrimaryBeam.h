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
