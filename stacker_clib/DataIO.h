#include <string>
#include <exception>
#include <stdexcept>

using std::string;
class Chunk;

#ifndef __DATA_H__
#define __DATA_H__

class fileException : public std::runtime_error {
	public:
		static const int OPEN = 0;
		static const int MOSAIC = 1;
		static const int HEADER_INFO_MISSING = 2;

	private:
		string message;
		int code;

	public:
		fileException(int code, string message) throw () : 
			std::runtime_error(message.c_str()), code(code){};
		~fileException() throw () {};

		int getCode() {return code; };
};

class DataIO
{
	public:
		virtual int nvis() = 0;
		virtual int readChunk(Chunk& chunk) = 0;
		virtual void writeChunk(Chunk& chunk) = 0;
		virtual int nPointings() = 0;
		virtual float xPhaseCentre(int fieldID) = 0;
		virtual float yPhaseCentre(int fieldID) = 0;
		virtual void setPhaseCentre(int fieldID, double x, double y) = 0;

};

#endif //inclusion guard

