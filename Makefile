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
	cp stacker_clib/*.so ~/.casa/ipython/stacker/stacker_clib
