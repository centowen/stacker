class Chunk;

#ifndef __DATA_H__
#define __DATA_H__

class DataIO
{
	public:
		virtual int nvis() = 0;
		virtual int readChunk(Chunk& chunk) = 0;
		virtual void writeChunk(Chunk& chunk) = 0;
		virtual int nPointings() = 0;
		virtual float xPhaseCentre(int id) = 0;
		virtual float yPhaseCentre(int id) = 0;
		virtual void setPhaseCentre(int id, double x, double y) = 0;

};

#endif //inclusion guard

