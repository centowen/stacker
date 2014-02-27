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
// A simple test program to read in uv coverage and calculate maximum and miniumum baseline.

#include <ms/MeasurementSets/MSTable.h>
#include <ms/MeasurementSets/MSColumns.h>
#include <ms/MeasurementSets/MeasurementSet.h>
#include <tables/Tables/TableDesc.h>
#include <casa/Arrays/VectorIter.h>
#include <casa/Arrays/MatrixIter.h>
#include <iostream>
#include <cmath>
#include <ctime>

using std::cout;
using std::endl;
using std::fabs;

using casa::TableDesc;
using casa::MeasurementSet;
using casa::ROMSColumns;
using casa::IPosition;

int main(int argc, char* argv[])
{
// 	TableDesc simpleDesc = MS::requiredTableDesc();

	MeasurementSet ms("PTG1.ms");
	MeasurementSet ms_out("tmp.ms");
	ROMSColumns msc(ms);
	MSColumns msc_out(ms_out);

// 	for(int i = 0; i < 1; i++)
// 		cout << "scanNumber ("<< i << "): " << msc.scanNumber()(i) << endl;

	double umax = -1.;
	double vmax = -1.;
	double wmax = -1.;
	cout << "nrow: " << msc.nrow() << endl;
	clock_t t1 = std::clock();
	casa::Array<double> a = msc.uvw().getColumn(); 
	IPosition start(1,0);
	IPosition d(1,1);

	cout << "shape: " << a.shape() << endl;
	casa::VectorIterator<double>* it = (casa::VectorIterator<double>*) a.makeIterator(1);
	casa::MatrixIterator<double>* datait = (casa::MatrixIterator<double>*) msc.data().getColumn().makeIterator(1);

	double mean = 0.;
// 	casa::Matrix<double> mean_mat(datait->matrix().shape());
	casa::Matrix<double> mean_mat(2,7);
 	for(int j = 0; j < mean_mat.ncolumn(); j++)
 		for(int i = 0; i < mean_mat.nrow(); i++)
 			mean_mat(i,j) = 0.;

	while(!it->pastEnd())
	{
		casa::Vector<double> b = it->vector();
		casa::Matrix<double> data = datait->matrix();

		double u = fabs(b[0]);
		double v = fabs(b[1]);
		double w = fabs(b[2]);

		mean += data(0,0);
		for(int j = 0; j < mean_mat.ncolumn(); j++)
			for(int i = 0; i < mean_mat.nrow(); i++)
				mean_mat(i,j) += data(i,j);

		if( umax < u) umax = u;
		if( vmax < v) vmax = v;
		if( wmax < w) wmax = w;

		it->next();
		datait->next();
	}
	mean /= msc.nrow();
	for(int j = 0; j < mean_mat.ncolumn(); j++)
		for(int i = 0; i < mean_mat.nrow(); i++)
			mean_mat(i,j) /= msc.nrow();

	clock_t t2 = std::clock();

	cout << t2 - t1 << " ticks " << endl;
	cout << double(*msc.uvw()(0)[0].data()) << endl;
	cout << "mean: " << mean << endl;
	cout << "mean_mat: " << mean_mat << endl;
	cout << "umax: " << umax << endl;
	cout << "vmax: " << vmax << endl;
	cout << "wmax: " << wmax << endl;

	for(int i = 0; i < 1; i++)
		cout << "antennas ("<< i << "): " << msc.antenna1()(i) << " - " << msc.antenna2()(i) << endl;

	cout << "datadescid ("<< 0 << "): " << msc.dataDescId()(0) << endl;
	}
