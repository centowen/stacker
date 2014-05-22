CASABASE=/home/lindroos/sw/casapy
CASACOREBASE=/home/lindroos/sw/casapy
#CASABASE=/astro/casa/casapy

# Use this for casa v < 4.2
# CASACORE_INCLUDE=-I$(CASABASE)/include/ -I$(CASABASE)/include/casacore/casa
# CASACORE_LIBPATH=$(CASABASE)/lib64
# CPPFLAGS=$(CASACORE_INCLUDE)
# LDFLAGS=-L$(CASACORE_LIBPATH) -lcasacore $(CASACORE_LIBPATH)/libgfortran.so.1 $(CASACORE_LIBPATH)/liblapack.so.3 $(CASACORE_LIBPATH)/libblas.so.3 -lwcs $(CASACORE_LIBPATH)/libcfitsio.so.0 -ldl 

# Use this for casa v. >= 4.2
CASACORE_INCLUDE=-I$(CASABASE)/include/ -I$(CASABASE)/casacore
CASACORE_LIBPATH=$(CASABASE)/lib64
CPPFLAGS=$(CASACORE_INCLUDE)
LDFLAGS=-L$(CASACORE_LIBPATH) -lcasa_ms -lcasa_images $(CASACORE_LIBPATH)/libgfortran.so.1 $(CASACORE_LIBPATH)/liblapack.so.3 $(CASACORE_LIBPATH)/libblas.so.3 -lwcs $(CASACORE_LIBPATH)/libcfitsio.so.0 -ldl 

default:
	g++ -static-libgcc -O3 -fPIC -c $(CPPFLAGS) $(DEBUG) -o stacker.o stacker.cpp
	g++ -static-libgcc -O3 -fPIC -c $(CPPFLAGS) -o Model.o Model.cpp
	g++ -static-libgcc -O3 -fPIC -c $(CPPFLAGS) -o MSComputer.o MSComputer.cpp
	g++ -static-libgcc -O3 -fPIC -c $(CPPFLAGS) -o StackChunkComputer.o StackChunkComputer.cpp
	g++ -static-libgcc -O3 -fPIC -c $(CPPFLAGS) -o ModsubChunkComputer.o ModsubChunkComputer.cpp
	g++ -static-libgcc -O3 -fPIC -c $(CPPFLAGS) -o PrimaryBeam.o PrimaryBeam.cpp
	g++ -static-libgcc -O3 -fPIC -c $(CPPFLAGS) -o Chunk.o Chunk.cpp
	g++ -static-libgcc -O3 -fPIC -c $(CPPFLAGS) -o msio.o msio.cpp
	g++ -static-libgcc -O3 -fPIC -c $(CPPFLAGS) -o Coords.o Coords.cpp
	g++ -static-libgcc -shared Coords.o Model.o msio.o Chunk.o PrimaryBeam.o ModsubChunkComputer.o StackChunkComputer.o MSComputer.o stacker.o -Wl,-soname,libstacker.so $(LDFLAGS) -o libstacker.so

clean:
	rm -f *.o *.last *.log *.log~ *.pyc

debug:
	DEBUG="-g" make default
