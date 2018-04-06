UNAME := $(shell uname)
ifeq ($(UNAME),Linux)
LIB_SUFFIX=so
endif
ifeq ($(UNAME),Darwin)
LIB_SUFFIX=dylib
endif

all:
	@echo
	@echo "cd stacker_clib and use scons to build"
	@echo "or \"make install\" to install prebuilt version for CASA"
	@echo

install:
	if [ -d ~/.casa/ipython/stacker ]; then rm -rf ~/.casa/ipython/stacker; fi
	mkdir -p ~/.casa/ipython/stacker
	cp -r  image __init__.py interval modsub pb uv ~/.casa/ipython/stacker
	mkdir -p ~/.casa/ipython/stacker/stacker_clib
	cp stacker_clib/*.$(LIB_SUFFIX) ~/.casa/ipython/stacker/stacker_clib
