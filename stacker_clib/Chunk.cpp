#include "Chunk.h"
#include <iostream>

Chunk::Chunk(int size)
{
	inVis = new Visibility[size];
	outVis = new Visibility[size];
	nvis = size;
}
Chunk::~Chunk()
{
	delete[] inVis, outVis;
	nvis = 0;
}

int Chunk::size()
{
	return nvis;
}

void Chunk::setSize(int size)
{
	if(size < 0)
		nvis = 0;
	else
	{
// 		delete[] inVis, outVis;
// 		inVis = new Visibility[size];
// 		outVis = new Visibility[size];
		nvis = size;
	}
}
