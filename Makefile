THIS_MAKEFILE_PATH:=$(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))
THIS_DIR:=$(shell cd $(dir $(THIS_MAKEFILE_PATH));pwd)

all:
	cd stacker_clib && make all

install:
	if [ ! -d ~/.casa/ipython ]; then mkdir -p ~/.casa/ipython; fi
	if [ -d ~/.casa/ipython/stacker ]; then rm -rf ~/.casa/ipython/stacker; fi
	cp -r $(THIS_DIR) ~/.casa/ipython
