CPPFLAGS= -DHAVE_KALLOC -std=c++11 -Wno-sign-compare -Wno-write-strings -Wno-unused-but-set-variable
LIBS= -lm -lz -lpthread

# Assuming git submodules were cloned previously
# If not, run "git submodule update --init --recursive" before "make"

all:winnowmap

winnowmap: MAKE_DIRS
	+$(MAKE) -C src
	$(CXX) $(CPPFLAGS) src/main.o -o bin/$@ -Lsrc -lminimap2 $(LIBS)
	+$(MAKE) -C src/meryl/src TARGET_DIR=$(shell pwd) 

MAKE_DIRS:
	@if [ ! -e bin ] ; then mkdir -p bin ; fi

clean:
	rm -rf bin
	rm -rf lib 
	+$(MAKE) clean -C src
	+$(MAKE) clean -C src/meryl/src

cleanw:
	rm -rf bin/winnowmap
	+$(MAKE) clean -C src
