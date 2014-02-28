#CASABASE=/home/lindroos/prog/casapy
CASABASE=/astro/casa/casapy
#CASABASE=/home/lindroos/prog/casapy-34.0.19988-002-64b
CASACORE_INCLUDE=-I$(CASABASE)/include/ -I$(CASABASE)/include/casacore/
CASACORE_LIBPATH=$(CASABASE)/lib64
#CASACORE_INCLUDE=-I/home/lindroos/prog/casapy-stable-33.0.16501-001-64b/include/ -I/home/lindroos/prog/casapy-stable-33.0.16501-001-64b/include/casacore/
#CASACORE_LIBPATH=-L/home/lindroos/prog/casapy-stable-33.0.16501-001-64b/lib64
ROOT_INCLUDE=
ROOT_LIBS=
#ROOT_INCLUDE=-I`root-config --incdir`
#ROOT_LIBS=`root-config --libs`


CPPFLAGS=$(CASACORE_INCLUDE) $(ROOT_INCLUDE)
#LDFLAGS=$(CASACORE_LIBPATH) $(ROOT_LIBS) -lcasacore -lgfortran.1 -llapack.3 -lblas.3 -lwcs -lcfitsio.0 -ldl 
LDFLAGS=-L$(CASACORE_LIBPATH) $(ROOT_LIBS) -lcasacore $(CASACORE_LIBPATH)/libgfortran.so.1 $(CASACORE_LIBPATH)/liblapack.so.3 $(CASACORE_LIBPATH)/libblas.so.3 -lwcs $(CASACORE_LIBPATH)/libcfitsio.so.0 -ldl 

default:
	g++ -O3 -fPIC -c $(CPPFLAGS) $(DEBUG) -o stacker.o stacker.cpp
	g++ -O3 -fPIC -c $(CPPFLAGS) -o Model.o Model.cpp
	g++ -O3 -fPIC -c $(CPPFLAGS) -o MSComputer.o MSComputer.cpp
	g++ -O3 -fPIC -c $(CPPFLAGS) -o StackChunkComputer.o StackChunkComputer.cpp
	g++ -O3 -fPIC -c $(CPPFLAGS) -o ModsubChunkComputer.o ModsubChunkComputer.cpp
	g++ -O3 -fPIC -c $(CPPFLAGS) -o PrimaryBeam.o PrimaryBeam.cpp
	g++ -O3 -fPIC -c $(CPPFLAGS) -o Chunk.o Chunk.cpp
	g++ -O3 -fPIC -c $(CPPFLAGS) -o msio.o msio.cpp
	g++ -O3 -fPIC -c $(CPPFLAGS) -o Coords.o Coords.cpp
	g++ -shared Coords.o Model.o msio.o Chunk.o PrimaryBeam.o ModsubChunkComputer.o StackChunkComputer.o MSComputer.o stacker.o -Wl,-soname,libstacker.so $(LDFLAGS) -o libstacker.so

clean:
	rm -f *.o *.last *.log *.log~ *.pyc

debug:
# 	g++ -fPIC -c $(CPPFLAGS) -o test.o stack.cpp
# 	g++ -O3 -fPIC -c $(CPPFLAGS) -o PrimaryBeam.o PrimaryBeam.cpp
# 	g++ PrimaryBeam.o test.o -Wl $(LDFLAGS) -o stack
	DEBUG="-g" make default
