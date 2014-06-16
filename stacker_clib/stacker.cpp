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

// includes/*{{{*/
#include "MSComputer.h"
#include "PrimaryBeam.h"
#include "MSPrimaryBeam.h"
#include "Model.h"
#include "Coords.h"
#include "ModsubChunkComputer.h"
#include "StackChunkComputer.h"
#include "definitions.h"
/*}}}*/

using std::cout;
using std::cerr;
using std::endl;


double cpp_stack(const char* infile, const char* outfile, 
			     int pbtype, char* pbfile, double* pbpar, int npbpar,
				 double* x, double* y, double* weight, int nstack);
void cpp_modsub(const char* infile, const char* outfile, 
	        const char* modelfile, const char* pbfile);

// Functions to interface with python module.
extern "C"{
	// Stacking function
	// Input arguments:
	// - infile: The input ms file.
	// - outfile: The output ms file, can be the same as input ms file.
	// - pbfile: A casa image of the primary beam, used to calculate primary beam correction.
	// - x: x coordinate of each stacking position (in radian).
	// - y: y coordinate of each stacking position (in radian).
	// - weight: weight of each stacking position.
	// - nstack: length of x, y and weight lists.
	// Returns average of all visibilities. Estimate of flux for point sources.
	//
	double stack(char* infile, char* outfile, 
			     int pbtype, char* pbfile, double* pbpar, int npbpar,
			     double* x, double* y, double* weight, int nstack)
	{
		double flux;
		flux = cpp_stack(infile, outfile, pbtype, pbfile, pbpar, npbpar, x, y, weight, nstack);
		return flux;
	};

	// Function to subtract model from uvdata
	// Input arguments:
	// - infile: The input ms file.
	// - outfile: The output ms file, can be the same as input ms file.
	// - infile: cl file with the model to be subtracted
	// - pbfile: A casa image of the primary beam, used to calculate primary beam correction.
	void modsub(char* infile, char* outfile, char* modelfile, char* pbfile)
	{
		cpp_modsub(infile, outfile, modelfile, pbfile);
	};
};

double cpp_stack(const char* infile, const char* outfile,
			     int pbtype, char* pbfile, double* pbpar, int npbpar,
				 double* x, double* y, double* weight, int nstack)/*{{{*/
{
	PrimaryBeam* pb;
	if(pbtype == PB_CONST)
		pb = (PrimaryBeam*)new ConstantPrimaryBeam;
	else if(pbtype == PB_MS)
		pb = (PrimaryBeam*)new MSPrimaryBeam(pbfile);
	else
		pb = (PrimaryBeam*)new ConstantPrimaryBeam;
	Coords coords(x, y, weight, nstack);
	StackChunkComputer* cc = new StackChunkComputer(&coords, pb);
	MSComputer* computer = new MSComputer((ChunkComputer*)cc, infile, outfile);

	float retval = computer->run();

	delete computer;
	delete cc;
	delete pb;

// 	Complex averageFlux =  float(1. / normSumVis) * sumVis;

	return retval;
}/*}}}*/

// Subtract a cl model from measurement set.
void cpp_modsub(const char* infile, const char* outfile, 
	        const char* modelfile, const char* pbfile) /*{{{*/
{
	PrimaryBeam* pb;// = new ImagePrimaryBeam(pbfile);
	pb = (PrimaryBeam*)new MSPrimaryBeam(pbfile);
	Model* model = new Model(modelfile);

	cout << "Pre making computer." << endl;

	ModsubChunkComputer* cc = new ModsubChunkComputer(model, pb);
	MSComputer* computer = new MSComputer(cc, infile, outfile);

	cout << "Pre running computer." << endl;

	computer->run();

	cout << "Post running computer." << endl;

	delete computer;
	delete cc;
	delete model;
	delete pb;
}
/*}}}*/
