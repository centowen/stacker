#include <iostream>
#include "PrimaryBeam.h"
#include "MSPrimaryBeam.h"
#include "Model.h"
#include "ModsubChunkComputer.h"
#include "definitions.h"

using std::cout;
using std::endl;

int main()
{
	PrimaryBeam* pb;// = new ImagePrimaryBeam(pbfile);
	pb = (PrimaryBeam*)new ConstantPrimaryBeam;
	cout << "Opening model." << endl;
	Model* model = new Model("../example/output/full.model");

	cout << "Creating cc." << endl;
	ModsubChunkComputer* cc = new ModsubChunkComputer(model, pb);
	int infiletype = FILE_TYPE_MS, outfiletype = FILE_TYPE_MS,
		infileoptions = 0, outfileoptions = 0;
// 	char infile[] = "../example/testdata.ms";
	char infile[] = "/home/lindroos/jobb/stacker/example/testdata.ms/";
	char outfile[] = "../example/output/residual.ms";
	cout << "Creating comp." << endl;
	MSComputer* computer = new MSComputer((ChunkComputer*)cc, 
			infiletype, infile, infileoptions,
			outfiletype, outfile, outfileoptions);
	cout << "Number of visibilities: " << computer->getMS()->nvis() << endl;
	cout << "Number of pointings: " << computer->getMS()->nPointings() << endl;

	return 0;
}

