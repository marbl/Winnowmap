CC= g++
CFLAGS= -g -Wall -O2 #-Wextra
CPPFLAGS= -DHAVE_KALLOC -std=c++11 -Wno-sign-compare -Wno-write-strings -Wno-unused-but-set-variable
PROG=		winnowmap
LIBS=		-lm -lz -lpthread

all:$(PROG)

winnowmap:
	+$(MAKE) -C src
	$(CC) $(CPPFLAGS) src/main.o -o $@ -Lsrc -lminimap2 $(LIBS)
	$(CC) $(CPPFLAGS) src/computeHighFreqKmers.cpp -o computeHighFreqKmers $(LIBS)

clean:
	rm -fr $(PROG) computeHighFreqKmers
	+$(MAKE) clean -C src

