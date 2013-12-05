// A simple test program to read in uv coverage and calculate maximum and miniumum baseline.


// includes/*{{{*/
#include <ms/MeasurementSets/MSTable.h>
#include <ms/MeasurementSets/MSColumns.h>
#include <ms/MeasurementSets/MeasurementSet.h>
#include <casa/Arrays/VectorIter.h>
#include <casa/Arrays/MatrixIter.h>
#include <images/Images/ImageOpener.h>
#include <images/Images/ImageInfo.h>
#include <lattices/Lattices/Lattice.h>
#include <lattices/Lattices/LatticeBase.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <tables/Tables/ExprNodeSet.h>
#include <tables/Tables/ExprNode.h>
#include <iostream>
#include <cmath>
#include <ctime>
#include <sstream>
#include <sys/stat.h>
#include <pthread.h>
#include <utility>
#include <stdio.h>
#include "ModsubChunkComputer.h"
#include "MSComputer.h"
#include "PrimaryBeam.h"
#include "Chunk.h"
#include "Model.h"
#include "time.h"
#include "definitions.h"/*}}}*/

using std::cout;
using std::cerr;
using std::endl;
using casa::C::pi;


void cpp_modsub(const char* msinfile, const char* msoutfile,
                  const char* modelfile, const char* pbfile);

// This is extern c to ensure no complications when called from outside library.
// Takes several input arguments:
//  - msinfile: The input ms file.
//  - msoutfile: The output ms file, can be the same as input ms file.
//  - modelfile: A path to a casa component list.
//  - pbfile: A casa image of the primary beam, used to calculate primary beam correction.
extern "C"{
	void modsub(char* msinfile, char* msoutfile, char* modelfile, char* pbfile)
	{
		cpp_modsub(msinfile, msoutfile, modelfile, pbfile);
	};
};

// Initiates actual stacking, called from the c version.
void cpp_modsub(const char* msinfile, const char* msoutfile, 
		        const char* modelfile, const char* pbfile) /*{{{*/
{
	PrimaryBeam pb(pbfile);


	Model model(modelfile);
	ModsubChunkComputer* cc = new ModsubChunkComputer(&model, &pb);
	MSComputer* computer = new MSComputer(cc, msinfile, msoutfile);

	computer->run();

	cout << "Time for cleanup (modsub)." << endl;

	delete computer;
	delete cc;

	cout << "Done with explicit cleanup, implicit to go. (modsub)." << endl;
}/*}}}*/


