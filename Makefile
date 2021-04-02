export CPPFLAGS= -g -Wall -O2 -DHAVE_KALLOC -fopenmp -std=c++11 -Wno-sign-compare -Wno-write-strings -Wno-unused-but-set-variable
export LIBS= -lm -lz -lpthread

all:winnowmap

winnowmap: MAKE_DIRS
	+$(MAKE) -e -C src
	$(CXX) $(CPPFLAGS)  src/main.o -o bin/$@ -Lsrc -lwinnowmap $(LIBS)

MAKE_DIRS:
	@if [ ! -e bin ] ; then mkdir -p bin ; fi

clean:
	rm -rf bin/winnowmap
	+$(MAKE) clean -C src
	rm -rf lib 
