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
#include "Coords.h"
#include <casa/Arrays/Array.h>
#include <ms/MeasurementSets/MeasurementSet.h>

using casa::MeasurementSet;
using casa::Array;

// Coords::Coords(const char* msfile, PrimaryBeam& pb, double* x, double* y, double* weight, int nstack)
Coords::Coords(double* x, double* y, double* weight, int nstack)
{
	x_raw = x;
	y_raw = y;
	weight_raw = weight;
	nStackPoints_raw = nstack;
	nPointings = 0;
	nStackPoints = NULL;
	omega_x = NULL;
	omega_y = NULL;
	omega_z = NULL;
	dx = NULL;
	dy = NULL;
	x = NULL;
	y = NULL;
	weight = NULL;
}


void Coords::computeCoords(DataIO* ms, PrimaryBeam& pb)
{
	nPointings = ms->nPointings();

    vector<float>* cx = new vector<float>[nPointings];
    vector<float>* cy = new vector<float>[nPointings];
    vector<float>* cweight = new vector<float>[nPointings];

    nStackPoints = new int[nPointings];
	nStackPointsVisible = 0;
    for(int i = 0; i < nPointings; i++)
        nStackPoints[i] = 0;

	for( int i = 0; i< nStackPoints_raw; i++ )
	{
		bool pointVisible = false;
		for(int fieldID = 0; fieldID < nPointings; fieldID++)
		{
            float dx = (float(x_raw[i]) - ms->xPhaseCentre(fieldID))*cos(y_raw[i]);
            float dy = asin(sin(float(y_raw[i]))/cos(dx)) - ms->yPhaseCentre(fieldID);
//          float dx = (x - x_phase_centre[fieldID])*cos(y);
//          float dy = asin(sin(y)) - y_phase_centre[fieldID];
			while(dx > 2*pi) dx -= 2*pi;
			while(dx < -2*pi) dx += 2*pi;

            if(pb.calc(dx,  dy )> 0.001)
            {
                cx[fieldID].push_back(float(x_raw[i]));
                cy[fieldID].push_back(float(y_raw[i]));
                if(weight_raw[i] <= tol)
					cweight[fieldID].push_back(1.);
				else
					cweight[fieldID].push_back(float(weight_raw[i]));
                nStackPoints[fieldID] ++;
				pointVisible = true;
            }
        }

		if(pointVisible)
			nStackPointsVisible ++;
    }

	for(int fieldID = 0; fieldID < nPointings; fieldID++)
		cout << "Number of stacking positions in field " << fieldID 
			 << ": " << nStackPoints[fieldID] << endl;

	omega_x = new float*[nPointings];
	omega_y = new float*[nPointings];
	omega_z = new float*[nPointings];
	dx = new float*[nPointings];
	dy = new float*[nPointings];
	this->x = new float*[nPointings];
	this->y = new float*[nPointings];
	this->weight = new float*[nPointings];

//     for(int fieldID = 0; fieldID < nPointings; fieldID++)
// 	{
//         cout << "Number of stacking positions in field " << fieldID << ": " << nStackPoints[fieldID] << endl;
// 	}

    for(int fieldID = 0; fieldID < nPointings; fieldID++)
    {
// 		cout << "field " << fieldID << ": " << endl;
        dx[fieldID] = new float[nStackPoints[fieldID]];
        dy[fieldID] = new float[nStackPoints[fieldID]];
        omega_x[fieldID] = new float[nStackPoints[fieldID]];
        omega_y[fieldID] = new float[nStackPoints[fieldID]];
        omega_z[fieldID] = new float[nStackPoints[fieldID]];
        this->x[fieldID] = new float[nStackPoints[fieldID]];
        this->y[fieldID] = new float[nStackPoints[fieldID]];
        this->weight[fieldID] = new float[nStackPoints[fieldID]];

        for(int i = 0; i < nStackPoints[fieldID]; i++)
        {
            this->x[fieldID][i] = cx[fieldID][i];
            this->y[fieldID][i] = cy[fieldID][i];
            this->weight[fieldID][i] = cweight[fieldID][i];

            dx[fieldID][i] = (this->x[fieldID][i] - ms->xPhaseCentre(fieldID))*cos(this->y[fieldID][i]);
            dy[fieldID][i] = asin(sin(this->y[fieldID][i])/cos(dx[fieldID][i])) - ms->yPhaseCentre(fieldID);
// 			std::cout << "pos " << i << ", dx, dy: " << dx << ", " << dy << endl;
			while(dx[fieldID][i] > 2*pi) dx[fieldID][i] -= 2*pi;
			while(dx[fieldID][i] < -2*pi) dx[fieldID][i] += 2*pi;
// 			cout << "(dx, dy) = (" << dx[fieldID][i] << ", " << dy[fieldID][i] << ")" << endl;
// 			std::cout << "pos " << i << ", dx, dy: " << dx << ", " << dy << endl;
// 
//          dx[fieldID][i] = (x[fieldID][i] - x_phase_centre[fieldID])*cos(y[fieldID][i]);
//          dy[fieldID][i] = asin(sin(y[fieldID][i])) - y_phase_centre[fieldID];

            omega_x[fieldID][i] = 2*pi*sin(dx[fieldID][i])/c;
            omega_y[fieldID][i] = 2*pi*sin(dy[fieldID][i])/c;
//             omega_z[fieldID][i] = 2*pi*(cos(sqrt(dx[fieldID][i]*dx[fieldID][i]+dy[fieldID][i]*dy[fieldID][i]))-1)/c;
            omega_z[fieldID][i] = 2*pi*(sqrt(1-dx[fieldID][i]*dx[fieldID][i]-dy[fieldID][i]*dy[fieldID][i])-1)/c;
// 			cout << "(omega_x, omega_y, omega_z) = (" 
// 				 << omega_x[fieldID][i] << ", " 
// 				 << omega_y[fieldID][i] << ", " 
// 				 << omega_z[fieldID][i] << ")" << endl;
        }

    }

	delete[] cx;
	delete[] cy;
	delete[] cweight;

	for(int i = 0; i < nStackPoints[0]; i++)
		cout << "coord " << i << ": " << omega_x[0][0] << ", " << omega_y[0][0] << "," << omega_z[0][0] << endl;
}

Coords::Coords(const char* coordfile)
{
    if(stat(coordfile, &statbuffer) == -1)
    {
        cerr << "Can not find coordinate file!" << endl;
        exit(-1);
    }
    ifstream file_csv(coordfile);

	vector<float> cx, cy, cweight;

	nPointings = 0;
	nStackPoints = NULL;
	omega_x = NULL;
	omega_y = NULL;
	omega_z = NULL;
	dx = NULL;
	dy = NULL;
	x = NULL;
	y = NULL;
	weight = NULL;

    string coordString("");
    getline ( file_csv, coordString);
    while ( !file_csv.eof() )
    {
        float x, y, weight;
        x = 0.;
        y = 0.;
        weight = 0.;

        stringstream coordReader(stringstream::in | stringstream::out);;
        coordReader.str("");

        coordString[coordString.find(",")] = ' ';
        coordString[coordString.find(",")] = ' ';
        coordReader << coordString;
        coordReader >> x;
        coordReader >> y;
        coordReader >> weight;
        coordReader.flush();

		cx.push_back(x);
		cy.push_back(y);
		if(weight <= tol) weight = 1.;
		cweight.push_back(weight);

        getline ( file_csv, coordString);
    }

	nStackPoints_raw = cx.size();
	x_raw = new double[nStackPoints_raw];
	y_raw = new double[nStackPoints_raw];
	weight_raw = new double[nStackPoints_raw];
	for(int i = 0; i < nStackPoints_raw; i++)
	{
		x_raw[i] = cx[i];
		y_raw[i] = cy[i];
		weight_raw[i] = cweight[i];
	}
}

Coords::~Coords()
{
	if(nPointings > 0)
	{
		for(int i = 0; i < nPointings; i++)
		{
			delete[] omega_x[i];
			delete[] omega_y[i];
			delete[] omega_z[i];
			delete[] dx[i];
			delete[] dy[i];
			delete[] x[i];
			delete[] y[i];
			delete[] weight[i];
		}
	}

	delete[] nStackPoints;
	delete[] omega_x;
	delete[] omega_y;
	delete[] omega_z;
	delete[] dx;
	delete[] dy;
	delete[] x;
	delete[] y;
	delete[] weight;
}
