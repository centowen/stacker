#include <casa/complex.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Vector.h>
#include <iostream>

#ifndef __CHUNK_H__
#define __CHUNK_H__

using casa::Complex;
using casa::Matrix;
using casa::Vector;

struct Visibility
{
	float u,v,w;
	double* freq;
	int fieldID, index, spw;

// 	float** data_real;
// 	float** data_imag;
// 	float* weight;

	Matrix<Complex> data;
	Vector<float> weight;

	public:
// 	Visibility()
// 	{
// 		data.resize(1,1);
// 		weight.resize(1);
// 	};
	~Visibility()
	{
	}
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
