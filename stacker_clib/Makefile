CASABASE=$(CASACORE_SW_DIR)# Path to CASA

#CASACORE_INCLUDE=-Iinclude
CASACORE_LIBPATH=$(CASABASE)/lib64
#CASACORE_INCLUDE=-I/home/lindroos/sw/casacore-1.7.0/
CASACORE_INCLUDE=-I/home/lindroos/sw/casacore-casapy4.2.2
#CASACORE_LIBPATH=/home/lindroos/sw/casacore-1.7.0/build/lib
CPPFLAGS=$(CASACORE_INCLUDE)
CASACORELIBS=-lcasa_scimath -lcasa_casa -lcasa_tables -lcasa_images -lcasa_ms -lcasa_components -lcasa_coordinates -lcasa_lattices -lcasa_measures -lcasa_lattices -lcasa_scimath_f -lcasa_fits -lcasa_msfits -lcasa_derivedmscal -lcasa_mirlib
# LDFLAGS=-L$(CASACORE_LIBPATH) -lcasa_ms -lcasa_images $(CASACORE_LIBPATH)/libgfortran.so.1 $(CASACORE_LIBPATH)/liblapack.so.3 $(CASACORE_LIBPATH)/libblas.so.3 -lwcs $(CASACORE_LIBPATH)/libcfitsio.so.0 -ldl 
LDFLAGS=-L$(CASACORE_LIBPATH) $(CASACORELIBS) $(CASACORE_LIBPATH)/libgfortran.so.1 $(CASACORE_LIBPATH)/liblapack.so.3 $(CASACORE_LIBPATH)/libblas.so.3 $(CASACORE_LIBPATH)/libwcs.so.4.7 $(CASACORE_LIBPATH)/libcfitsio.so.0 -ldl
#LDFLAGS=-L$(CASACORE_LIBPATH) $(CASACORELIBS) -lgfortran -llapack -lblas -lwcs -lcfitsio -ldl
OBJECTS=Coords.o Model.o DataIOFits.o msio.o Chunk.o PrimaryBeam.o MSPrimaryBeam.o ModsubChunkComputer.o StackChunkComputer.o MSComputer.o stacker.o
TESTOBJECTS=Coords.o Model.o DataIOFits.o msio.o Chunk.o PrimaryBeam.o MSPrimaryBeam.o ModsubChunkComputer.o StackChunkComputer.o MSComputer.o stacker_test.o
DEPS = $(wildcard *.h)

all: $(OBJECTS)
	g++ -static-libgcc -shared $(OBJECTS) -Wl,-soname,libstacker.so $(LDFLAGS) -o libstacker.so

test: $(TESTOBJECTS)
	g++ $(TESTOBJECTS) $(LDFLAGS) $(LDFLAGS) -o stacker_test

%.o: %.cpp ${DEPS}
	g++ -fPIC -c $(CPPFLAGS) $(DEBUG) -o $@ $<

clean:
	rm -f *.o libstacker.so *.last *.log *.log~ *.pyc

debug:
	DEBUG="-g" make all
	DEBUG="-g" make test
