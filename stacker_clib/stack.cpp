// includes/*{{{*/
#include <ms/MeasurementSets/MSTable.h>
#include <ms/MeasurementSets/MSColumns.h>
#include <ms/MeasurementSets/MeasurementSet.h>
#include <tables/Tables/TableDesc.h>
#include <casa/Arrays/VectorIter.h>
#include <casa/Arrays/MatrixIter.h>
#include <images/Images/ImageOpener.h>
#include <images/Images/ImageInfo.h>
#include <images/Images/ImageInterface.h>
#include <lattices/Lattices/Lattice.h>
#include <lattices/Lattices/LatticeBase.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <tables/Tables/ExprNodeSet.h>
#include <tables/Tables/ExprNode.h>
#include <casa/OS/Directory.h>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <pthread.h>
#include <queue>
#include <utility>
#include <string>
#include <string.h>
#include <stdio.h>
#include "StackChunkComputer.h"
#include "MSComputer.h"
#include "PrimaryBeam.h"
#include "Chunk.h"
#include "Coords.h"
#include "time.h"
#include "definitions.h"/*}}}*/

// using statements/*{{{*/
using std::cout;
using std::cerr;
using std::endl;
using std::fabs;
using std::exp;
using std::fstream;
using std::ifstream;
using std::getline;
using std::stringstream;
using std::vector;
using std::cos;
using std::sin;
using std::asin;
using std::queue;
using std::pair;
using std::string;

using casa::C::pi;

using casa::TableDesc;
using casa::IPosition;
using casa::Complex;
using casa::Float;
using casa::ImageInterface;
using casa::ImageOpener;
using casa::CoordinateSystem;
using casa::Directory;/*}}}*/

string to_string(int x)
{
	return dynamic_cast< std::ostringstream & >( \
		        ( std::ostringstream() << std::dec << x ) ).str();
};

void calculate_chunk_average(Chunk& chunk);

double cpp_stack(const char* msinfile, const char* msoutfile, const char* pbfile, double* x, double* y, double* weight, int nstack);

// This is extern c to ensure no complications when called from outside library.
// Takes several input arguments:
//  - msinfile: The input ms file.
//  - msoutfile: The output ms file, can be the same as input ms file.
//  - pbfile: A casa image of the primary beam, used to calculate primary beam correction.
//  - x: x coordinate of each stacking position.
//  - y: y coordinate of each stacking position.
//  - weight: weight of each stacking position.
//  - nstack: length of x, y and weight lists.
//  - inStackingMode: Allows multiply each visibility with a psf function.
//  - inBeam: Parameter for inStackingMode psf function.
//  Returns average of all visibilities. Only useful is the case of simple data.
extern "C"{
	double stack(char* msinfile, char* msoutfile, char* pbfile, double* x, 
			     double* y, double* weight, int nstack, int inStackingMode, 
				 double inBeam)
	{
		double flux;

// 		stackingMode = inStackingMode;
// 		if(stackingMode)
// 			beam = inBeam;
// 		else
// 			beam = 0.;

		flux = cpp_stack(msinfile, msoutfile, pbfile, x, y, weight, nstack);

		return flux;
	};
};

// Initiates actual stacking, called from the c version.
double cpp_stack(const char* msinfile, const char* msoutfile, const char* pbfile,
		double* x, double* y, double* weight, int nstack)/*{{{*/
{
	cout << "In c++!" << endl;
	PrimaryBeam pb(pbfile);
	Coords coords(x, y, weight, nstack);
	StackChunkComputer* cc = new StackChunkComputer(&coords, &pb);
	cout << "made cc" << endl;
	MSComputer* computer = new MSComputer(cc, msinfile, msoutfile);
	cout << "made main computer" << endl;

	computer->run();
	cout << "done running main computer" << endl;

// 	msdealer.setPhaseCenter(0, 0., 0.);

	delete computer;
	delete cc;

// 	Complex averageFlux =  float(1. / normSumVis) * sumVis;

	return 0.0;
}/*}}}*/



// Computes a running average to return at the end.
// Quick estimate for stacked flux as long as the source in a point source.
// void calculate_chunk_average(Chunk& chunk)
// {
// 	Complex sumVisLocal(float(0.),float(0.));
// 	float normSumVisLocal = 0.;
// 	for( int i = 0; i < chunk.size(); i++)
// 	{
// 		Complex aVis(0.,0.);
// 		for( int row = 0; row < chunk.outVis[i].data.nrow(); row ++)
// 			for( int column = 0; column < chunk.outVis[i].data.ncolumn(); column ++)
// 			{
// 				aVis += chunk.outVis[i].data(row, column);
// 			}
// 		aVis = float(1./(chunk.outVis[i].data.nrow()*chunk.outVis[i].data.ncolumn()))*aVis;
// 		sumVisLocal += float(chunk.outVis[i].weight(0))*aVis;
// 		normSumVisLocal += float(chunk.outVis[i].weight(0));
// 	}
// 	pthread_mutex_lock(&mutex);
// 	sumVis += sumVisLocal ;
// 	normSumVis += normSumVisLocal ;
// 	pthread_mutex_unlock(&mutex);
// }
// 
// // Should indicate a psf to be used for stacking. Needs to be specfied in the
// // uv plane. The simplest is of course dd = dd => real point. But this can 
// // also be used to take things such as point in middle minus a circle around.
// Complex psf(double vis_real, double vis_imag, double u, double v, double w, double dx, double dy, double freq)
// {
// 	// Phi = -2*pi*freq/c*(u*dx+v*dy)
// 	// or more correct
// 	// Phi = -2*pi*freq/c*(u*sin(dx)+v*sin(dy)+(cos(sqrt(dx**2+dy**2))-1)*w)
// 	double Phi_u, Phi_v, Phi, Phi_v_h, Phi_u_h;
// 
// 	if(stackingMode == 1)
// 	{
// 		Phi_u = -2*pi*freq/casa::C::c*u*beam/2.;
// 		Phi_v = -2*pi*freq/casa::C::c*v*beam/2.;
// 		Phi_u_h = -2*pi*freq/casa::C::c*(u+v)/sqrt(2)*beam/2.;
// 		Phi_v_h = -2*pi*freq/casa::C::c*(u-v)/sqrt(2)*beam/2.;
// 		vis_real = (1-(0.25*cos(Phi_u)+0.25*cos(Phi_v)+0.25*cos(Phi_u_h)+0.25*cos(Phi_v_h)))*vis_real;
// 	}
// 	else if(stackingMode == 2)
// 	{
// 		Phi = -2*pi*freq/casa::C::c/2.*(u*beam*cos(dx) + sin(dx)*beam*w);
// 		vis_real = (cos(Phi))*vis_real;
// 		vis_imag = (sin(Phi))*vis_imag;
// // 		vis_real = (1-cos(Phi))*vis_real;
// // 		vis_imag = (1-sin(Phi))*vis_real;
// 	}
// 	else if(stackingMode == 3)
// 	{
// 		Phi = 2*pi*freq/casa::C::c*beam*sqrt(u*u+v*v);
// 		vis_real = j0(Phi)*vis_real;
// 	}
// 
// 	return Complex(vis_real, vis_imag);
// }
