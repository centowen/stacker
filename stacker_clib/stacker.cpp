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
#include "config.h"
#ifdef USE_CUDA
#include "StackChunkComputerGpu.h"
#include "ModsubChunkComputerGpu.h"
#include "StackMCCCGpu.h"
#endif
/*}}}*/

using std::cout;
using std::cerr;
using std::endl;


void cpp_stack_mc(int infiletype, const char* infile, int infileoptions, 
                  int pbtype, char* pbfile, double* pbpar, int npbpar,
                  double* x, double* y, double* weight, int nstack, int nmc,
                  char** modelfiles, 
                  double* res_flux, double* res_weight, double* limits, int nbin,
                  bool use_cuda = true);
double cpp_stack(int infiletype, const char* infile, int infileoptions, 
                 int outfiletype, const char* outfile, int outfileoptions, 
                 int pbtype, char* pbfile, double* pbpar, int npbpar,
                 double* x, double* y, double* weight, int nstack,
                 bool use_cuda = false);
void cpp_modsub(int infiletype, const char* infile, int infileoptions, 
                int outfiletype, const char* outfile, int outfileoptions, 
                const char* modelfile,
                int pbtype, const char* pbfile, double* pbpar, int npbpar,
				bool subtract = true, bool use_cuda = false,
				const bool selectField=false, const char* field="");

// Functions to interface with python module.
extern "C"{/*{{{*/
	// Stacking function/*{{{*/
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
	             double* x, double* y, double* weight, int nstack,
	             bool use_cuda = false)
	{
		double flux;
		flux = cpp_stack(infiletype, infile, infileoptions, 
		                 outfiletype, outfile, outfileoptions,
		                 pbtype, pbfile, pbpar, npbpar, 
		                 x, y, weight, nstack, use_cuda);
		return flux;
	};/*}}}*/

	// Stacking Monte-Carlo noise esimate function/*{{{*/
	// Input arguments:
	// - infile: The input ms file.
	// - pbparameters: Generate using the appropriate function in the pb module.
	// - x: x coordinate of each stacking position (in radian).
	// - y: y coordinate of each stacking position (in radian).
	// - weight: weight of each stacking position.
	// - nstack: length of x, y and weight lists.
	// - nmc: Number of Monte-Carlo iterations.
	// - modelfiles: Models to add for each Monte-Carlo iteration.
	// - res_flux: Array to write result to (must be nbin*nmc long).
	// - res_weight: Array to write result to (must be nbin*nmc long).
	// - nbin: Number of bins to calculate flux in.
	// - use_cuda: Switch to 
	// Returns average of all visibilities. Estimate of flux for point sources.
	//
	void stack_mc(int infiletype, const char* infile, int infileoptions, 
	              int pbtype, char* pbfile, double* pbpar, int npbpar,
	              double* x, double* y, double* weight, int nstack, int nmc,
	              char** modelfiles, 
	              double* res_flux, double* res_weight, double* bins, int nbin,
	              bool use_cuda = true)
	{
		cpp_stack_mc(infiletype, infile, infileoptions,
		             pbtype, pbfile, pbpar, npbpar,
		             x, y, weight, nstack, nmc,
		             modelfiles,
		             res_flux, res_weight, bins, nbin,
		             use_cuda);
// 		int i;
// 		for(i = 0; i < nmc; i++)
// 			cout << modelfiles[i] << endl;
	};/*}}}*/

	// Function to subtract model from uvdata/*{{{*/
	// Input arguments:
	// - infile: The input ms file.
	// - outfile: The output ms file, can be the same as input ms file.
	// - infile: cl file with the model to be subtracted
	// - pbfile: A casa image of the primary beam, used to calculate primary beam correction.
	void modsub(int infiletype, char* infile, int infileoptions, 
	            int outfiletype, char* outfile, int outfileoptions,
	            char* modelfile, 
	            int pbtype, const char* pbfile, double* pbpar, int npbpar,
	            bool subtract = true, bool use_cuda = false,
				const bool selectField = false, const char* field = "")
	{
		cpp_modsub(infiletype, infile, infileoptions, 
		           outfiletype, outfile, outfileoptions,
		           modelfile, 
		           pbtype, pbfile, pbpar, npbpar,
		           subtract, use_cuda, selectField, field);
	};/*}}}*/
};/*}}}*/

void cpp_stack_mc(int infiletype, const char* infile, int infileoptions, /*{{{*/
                  int pbtype, char* pbfile, double* pbpar, int npbpar,
                  double* x, double* y, double* weight, int nstack, int nmc,
                  char** modelfiles, 
                  double* res_flux, double* res_weight, double* bins, int nbin,
                  bool use_cuda)
{
#ifdef USE_CUDA
	PrimaryBeam* pb;
	if(pbtype == PB_CONST)
		pb = (PrimaryBeam*)new ConstantPrimaryBeam;
	else if(pbtype == PB_MS)
	{
		pb = (PrimaryBeam*)new MSPrimaryBeam(pbfile);
	}
	else
		pb = (PrimaryBeam*)new ConstantPrimaryBeam;

	Coords** coordlists = new Coords*[nmc];
	Model** models = new Model*[nmc];

	cout << "bins: ";
	for(int i = 0; i < nbin+1; i++)
	{
		cout << bins[i] << ", ";
	}
	cout << endl;

	for(int i = 0; i < nmc; i++)
	{
		coordlists[i] = new Coords(&x[i*nstack], &y[i*nstack],
				                   &weight[i*nstack], nstack);
		if(strcmp(modelfiles[i], "") != 0)
			models[i] = new Model(modelfiles[i], false);
		else
			models[i] = NULL;
	}

	StackMCCCGpu* cc;
	int n_thread = N_THREAD;
	if(use_cuda)
	{
		cout << "Creating cc." << endl;
		cc = new StackMCCCGpu(coordlists, models, pb, nmc, bins, nbin);
		n_thread = 1;
	}
	else
	{
		cout << "Non CUDA support is not compiled. Recompile to enable." 
		     << endl;
		return;
	}

	MSComputer* computer;
	try
	{
		cout << "Creating computer." << endl;
		cout << "infile: " << infile << " with type " << infiletype << endl;
		computer = new MSComputer((ChunkComputer*)cc, 
								  infiletype, infile, infileoptions,
								  FILE_TYPE_NONE, "", 0,
								  n_thread);
		cout << "Computer created." << endl;
		computer->run();
	}
	catch(fileException e)
	{
		std::cerr << e.what() << std::endl;
	}
// 	double averageFlux =  (double)cc->flux();
	cout << "Time to copy results." << endl;
	double* bin_flux = cc->get_flux();
	double* bin_weight = cc->get_weight();
	for(int i = 0; i < nmc*nbin; i++)
	{
		res_flux[i] = bin_flux[i];
		res_weight[i] = bin_weight[i];
	}
	delete[] bin_flux;
	delete[] bin_weight;

	cout << "time to delete stuff" << endl;
	for(int i = 0; i < nmc; i++)
	{
		delete coordlists[i];
		if(models[i] != NULL)
			delete models[i];
	}
	delete[] coordlists;
	delete[] models;
	delete computer;
	delete cc;
	delete pb;
#else
	cout << "CUDA support is not compiled. Recompile to enable." << endl;
	return;
#endif
}/*}}}*/
double cpp_stack(int infiletype, const char* infile, int infileoptions, /*{{{*/
                 int outfiletype, const char* outfile, int outfileoptions, 
			     int pbtype, char* pbfile, double* pbpar, int npbpar,
				 double* x, double* y, double* weight, int nstack,
				 bool use_cuda)
{
	PrimaryBeam* pb;
	if(pbtype == PB_CONST)
		pb = (PrimaryBeam*)new ConstantPrimaryBeam;
	else if(pbtype == PB_MS)
	{
		pb = (PrimaryBeam*)new MSPrimaryBeam(pbfile);
	}
	else
		pb = (PrimaryBeam*)new ConstantPrimaryBeam;

	Coords coords(x, y, weight, nstack);
	ChunkComputer* cc;
	int n_thread = N_THREAD;
	if(use_cuda)
	{
#ifdef USE_CUDA
		cc = (ChunkComputer*) new StackChunkComputerGpu(&coords, pb);
		n_thread = 1;
#else
		cout << "CUDA support is not compiled. Recompile to enable." << endl;
		return 0.;
#endif
	}
	else
	{
		cc = (ChunkComputer*) new StackChunkComputer(&coords, pb);
	}

	MSComputer* computer;
	try
	{
		computer = new MSComputer(cc, 
								  infiletype, infile, infileoptions,
								  outfiletype, outfile, outfileoptions,
								  n_thread);
		computer->run();
	}
	catch(fileException e)
	{
		std::cerr << e.what() << std::endl;
	}
// 	double averageFlux =  (double)cc->flux();

	delete computer;
	delete cc;
	delete pb;


// 	return averageFlux;
	return 0.;
}/*}}}*/

// Subtract a cl model from measurement set.
void cpp_modsub(int infiletype, const char* infile, int infileoptions, /*{{{*/
                int outfiletype, const char* outfile, int outfileoptions, 
                const char* modelfile, 
                int pbtype, const char* pbfile, double* pbpar, int npbpar,
                bool subtract, bool use_cuda,
				const bool selectField, const char* field)
{
	PrimaryBeam* pb;// = new ImagePrimaryBeam(pbfile);
	if(pbtype == PB_CONST)
		pb = (PrimaryBeam*)new ConstantPrimaryBeam;
	else if(pbtype == PB_MS)
		pb = (PrimaryBeam*)new MSPrimaryBeam(pbfile);
	else
		pb = (PrimaryBeam*)new ConstantPrimaryBeam;

	Model* model = new Model(modelfile, subtract);


	cout << "subtract = " << subtract << endl;
	ChunkComputer* cc;
	int n_thread = N_THREAD;
	if(use_cuda)
	{
#ifdef USE_CUDA
		cout << "Going to use cuda for computations." << endl;
		cc = (ChunkComputer*) new ModsubChunkComputerGpu(model, pb);
		n_thread = 1;
#else
		cout << "CUDA support is not compiled. Recompile to enable." << endl;
		return;
#endif
	}
	else
	{
		cc = (ChunkComputer*) new ModsubChunkComputer(model, pb);
	}
	MSComputer* computer = NULL;
	try
	{
		computer = new MSComputer(cc, 
		                          infiletype, infile, infileoptions,
		                          outfiletype, outfile, outfileoptions,
		                          n_thread, selectField, field);
		computer->run();

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
