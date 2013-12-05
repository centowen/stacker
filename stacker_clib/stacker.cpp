// Library to stack and modsub ms data.

// includes/*{{{*/
#include "MSComputer.h"
#include "PrimaryBeam.h"
#include "Model.h"
#include "Coords.h"
#include "ModsubChunkComputer.h"
#include "StackChunkComputer.h"
/*}}}*/

using std::cout;
using std::cerr;
using std::endl;


double cpp_stack(const char* msinfile, const char* msoutfile, const char* pbfile,
		double* x, double* y, double* weight, int nstack);
void cpp_modsub(const char* msinfile, const char* msoutfile, 
	        const char* modelfile, const char* pbfile);

extern "C"{
	double stack(char* msinfile, char* msoutfile, char* pbfile, double* x, 
			     double* y, double* weight, int nstack)
	{
		double flux;
		flux = cpp_stack(msinfile, msoutfile, pbfile, x, y, weight, nstack);
		return flux;
	};

	void modsub(char* msinfile, char* msoutfile, char* modelfile, char* pbfile)
	{
		cpp_modsub(msinfile, msoutfile, modelfile, pbfile);
	};
};

double cpp_stack(const char* msinfile, const char* msoutfile, const char* pbfile,
		double* x, double* y, double* weight, int nstack)/*{{{*/
{
	PrimaryBeam* pb = new PrimaryBeam(pbfile);
	Coords coords(x, y, weight, nstack);
	StackChunkComputer* cc = new StackChunkComputer(&coords, pb);
	MSComputer* computer = new MSComputer(cc, msinfile, msoutfile);

	computer->run();

	delete computer;
	delete cc;
	delete pb;

// 	Complex averageFlux =  float(1. / normSumVis) * sumVis;

	return 0.0;
}/*}}}*/

// Subtract a cl model from measurement set.
void cpp_modsub(const char* msinfile, const char* msoutfile, 
	        const char* modelfile, const char* pbfile) /*{{{*/
{
	PrimaryBeam* pb = new PrimaryBeam(pbfile);
	Model* model = new Model(modelfile);

	cout << "Pre making computer." << endl;

	ModsubChunkComputer* cc = new ModsubChunkComputer(model, pb);
	MSComputer* computer = new MSComputer(cc, msinfile, msoutfile);

	cout << "Pre running computer." << endl;

	computer->run();

	cout << "Post running computer." << endl;

	delete computer;
	delete cc;
	delete model;
	delete pb;
}
/*}}}*/
