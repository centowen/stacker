// #include <casa/complex.h>
// #include <casa/Arrays/Matrix.h>
// #include <casa/Arrays/Vector.h>
#include <iostream>

#ifndef __CHUNK_H__
#define __CHUNK_H__

// using casa::Complex;
// using casa::Matrix;
// using casa::Vector;

struct Visibility
{
	float u,v,w;
	double* freq;
	float *data_real, *data_imag, *weight;
	int nstokes, nchan;
	int fieldID, index, spw;

public:
	Visibility();
	~Visibility();
};

class Chunk
{
private:
	int nvis;

public:
	Visibility *inVis, *outVis;
	Chunk(int size);
	~Chunk();
	int size();
	void setSize(int size);
};

#endif // end of inclusion guard
