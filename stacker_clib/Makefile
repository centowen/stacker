CASABASE=$(CASACORE_SW_DIR) # Path to CASA
CASACOREBASE=$(CASACORE_SW_DIR)/casacore # Path to CASACORE

CASACORE_INCLUDE=-I$(CASABASE)/include/ -I$(CASACOREBASE)
CASACORE_LIBPATH=$(CASABASE)/lib64
CPPFLAGS=$(CASACORE_INCLUDE)
LDFLAGS=-L$(CASACORE_LIBPATH) -lcasa_ms -lcasa_images $(CASACORE_LIBPATH)/libgfortran.so.1 $(CASACORE_LIBPATH)/liblapack.so.3 $(CASACORE_LIBPATH)/libblas.so.3 -lwcs $(CASACORE_LIBPATH)/libcfitsio.so.0 -ldl 
OBJECTS=Coords.o Model.o DataIOFits.o msio.o Chunk.o PrimaryBeam.o MSPrimaryBeam.o ModsubChunkComputer.o StackChunkComputer.o MSComputer.o stacker.o
DEPS = $(wildcard *.h)

all: $(OBJECTS)
	g++ -static-libgcc -shared $(OBJECTS) -Wl,-soname,libstacker.so $(LDFLAGS) -o libstacker.so

%.o: %.cpp ${DEPS}
	g++ -static-libgcc -O3 -fPIC -c $(CPPFLAGS) $(DEBUG) -o $@ $<

clean:
	rm -f *.o *.last *.log *.log~ *.pyc

debug:
	DEBUG="-g" make all
