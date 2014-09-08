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
#include "DataIO.h"
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


double cpp_stack(int infiletype, const char* infile, int infileoptions, 
                 int outfiletype, const char* outfile, int outfileoptions, 
			     int pbtype, char* pbfile, double* pbpar, int npbpar,
				 double* x, double* y, double* weight, int nstack);
void cpp_modsub(int infiletype, const char* infile, int infileoptions, 
                int outfiletype, const char* outfile, int outfileoptions, 
				const char* modelfile,
				int pbtype, const char* pbfile, double* pbpar, int npbpar);

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
	double stack(int infiletype, const char* infile, int infileoptions, 
                 int outfiletype, const char* outfile, int outfileoptions, 
			     int pbtype, char* pbfile, double* pbpar, int npbpar,
			     double* x, double* y, double* weight, int nstack)
	{
		double flux;
		flux = cpp_stack(infiletype, infile, infileoptions, 
				         outfiletype, outfile, outfileoptions,
						 pbtype, pbfile, pbpar, npbpar, 
						 x, y, weight, nstack);
		return flux;
	};

	// Function to subtract model from uvdata
	// Input arguments:
	// - infile: The input ms file.
	// - outfile: The output ms file, can be the same as input ms file.
	// - infile: cl file with the model to be subtracted
	// - pbfile: A casa image of the primary beam, used to calculate primary beam correction.
	void modsub(int infiletype, char* infile, int infileoptions, 
			    int outfiletype, char* outfile, int outfileoptions,
				char* modelfile, 
				int pbtype, const char* pbfile, double* pbpar, int npbpar)
	{
		cpp_modsub(infiletype, infile, infileoptions, 
				   outfiletype, outfile, outfileoptions,
				   modelfile, 
				   pbtype, pbfile, pbpar, npbpar);
	};
};

double cpp_stack(int infiletype, const char* infile, int infileoptions, 
                 int outfiletype, const char* outfile, int outfileoptions, 
			     int pbtype, char* pbfile, double* pbpar, int npbpar,
				 double* x, double* y, double* weight, int nstack)/*{{{*/
{
	PrimaryBeam* pb;
	if(pbtype == PB_CONST)
		pb = (PrimaryBeam*)new ConstantPrimaryBeam;
	else if(pbtype == PB_MS)
	{
		pb = (PrimaryBeam*)new MSPrimaryBeam(pbfile);
		cout << "Using ms pb model." << endl;
	}
	else
		pb = (PrimaryBeam*)new ConstantPrimaryBeam;
	cout << "nstack = " << nstack << endl;
	Coords coords(x, y, weight, nstack);
	StackChunkComputer* cc = new StackChunkComputer(&coords, pb);
	MSComputer* computer = new MSComputer((ChunkComputer*)cc, 
			                              infiletype, infile, infileoptions,
										  outfiletype, outfile, outfileoptions);

	float retval = computer->run();
	double averageFlux =  (double)cc->flux();

	delete computer;
	delete cc;
	delete pb;


	return averageFlux;
}/*}}}*/

// Subtract a cl model from measurement set.
void cpp_modsub(int infiletype, const char* infile, int infileoptions, 
                int outfiletype, const char* outfile, int outfileoptions, 
				const char* modelfile, 
				int pbtype, const char* pbfile, double* pbpar, int npbpar) /*{{{*/
{
	PrimaryBeam* pb;// = new ImagePrimaryBeam(pbfile);
	if(pbtype == PB_CONST)
		pb = (PrimaryBeam*)new ConstantPrimaryBeam;
	else if(pbtype == PB_MS)
		pb = (PrimaryBeam*)new MSPrimaryBeam(pbfile);
	else
		pb = (PrimaryBeam*)new ConstantPrimaryBeam;

	Model* model = new Model(modelfile);

	cout << "Pre making computer." << endl;

	ModsubChunkComputer* cc = new ModsubChunkComputer(model, pb);
	MSComputer* computer = NULL;
	try
	{
		computer = new MSComputer((ChunkComputer*)cc, 
											infiletype, infile, infileoptions,
											outfiletype, outfile, outfileoptions);
		cout << "Pre running computer." << endl;

		computer->run();

		cout << "Post running computer." << endl;
	}
	catch(fileException e)
	{
		std::cerr << e.what() << std::endl;
	}


	delete computer;
	delete cc;
	delete model;
	delete pb;
}
/*}}}*/
