#include "types.H"
#include "files.H"
#include "system.H"
#include "math.H"

#include <vector>
#include <deque>

class cpBuf {
public:
  cpBuf(uint64 maxsize) {
    _bLen = 0;
    _bMax = maxsize;
    _b    = new uint8 [_bMax];
  }

  ~cpBuf() {
    delete [] _b;
  }

  uint64  read(FILE *F) {
    return _bLen = merylutil::loadFromFile(_b, "stuff", _bMax, F, false);
  }


  uint64  _bMax = 0;
  uint64  _bLen = 0;
  uint8  *_b    = nullptr;
};



class catState {
public:
  catState() {
  }
  ~catState() {
  };

  cpBuf      *read(cpBuf *b = nullptr) {

    if ((inPath == nullptr) && (inFiles.size() == 0)) {
      delete b;
      return nullptr;
    }

    if (inPath == nullptr) {
      inFile = merylutil::openInputFile(inPath = inFiles.front());
      inFiles.pop_front();
    }

    if (b == nullptr)
      b = new cpBuf(blockSize);

    if (b->read(inFile) == 0) {
      merylutil::closeFile(inFile);

      inPath = nullptr;
      inFile = nullptr;

      return read(b);
    }

    return b;
  };

  void        write(cpBuf *b) {
    merylutil::writeToFile(b->_b, "stuff", b->_bLen, stdout);

    delete b;
  }

  //  Command line parameters.
  std::deque<char const *>   inFiles;

  uint64      maxMemory   = 0;
  uint64      blockSize   = 0;
  uint64      queueLength = 0;

  //  State;
  char const *inPath = nullptr;
  FILE       *inFile = nullptr;
};




void *bufReader(void *G) {
  catState *g = (catState *)G;
  return g->read();
}

void  bufWorker(void *G, void *T, void *S) {
}

void  bufWriter(void *G, void *S) {
  catState *g = (catState *)G;
  cpBuf    *s = (cpBuf    *)S;
  g->write(s);
}

void  bufStatus(void *G, uint64 numberLoaded, uint64 numberComputed, uint64 numberOutput) {
  catState *g = (catState *)G;
}




int
main(int argc, char **argv) {
  catState                   *g  = new catState();
  sweatShop                  *ss = new sweatShop(bufReader, bufWorker, bufWriter, bufStatus);
  std::vector<char const *>   e;

  for (int arg=1; arg<argc; arg++) {
    if      (strcmp(argv[arg], "-") == 0)
      g->inFiles.push_back(argv[arg]);

    else if (strcmp(argv[arg], "-queuelength") == 0)
      decodeInteger(argv[++arg], 0, 0, g->queueLength, e);

    else if (strcmp(argv[arg], "-blocksize") == 0)
      decodeInteger(argv[++arg], 0, 0, g->blockSize, e);

    else if (strcmp(argv[arg], "-maxmemory") == 0)
      decodeInteger(argv[++arg], 0, 0, g->maxMemory, e);

    else if (merylutil::fileExists(argv[arg]))
      g->inFiles.push_back(argv[arg]);

    else
      sprintf(e, "ERROR: '%s' is neither an input file nor stdin ('-').\n", argv[arg]);
  }

  bool  m = (g->maxMemory   != 0);   bool M = !m;   uint64 defm = 1llu << 30;
  bool  b = (g->blockSize   != 0);   bool B = !b;   uint64 defb = 1llu << 16;
  bool  q = (g->queueLength != 0);   bool Q = !q;   uint64 defq = 1llu << 14;

  if      (M && B && Q) { g->blockSize = defb;  g->queueLength = defq; }
  else if (M && B && q) { g->blockSize = defb;                         }
  else if (M && b && Q) {                       g->queueLength = defq; }
  else if (M && b && q) { /* perfect! */                               }

  else if (m && B && Q) { g->blockSize = defb;
                          g->queueLength = g->maxMemory / g->blockSize;   }
  else if (m && B && q) { g->blockSize   = g->maxMemory / g->queueLength; }
  else if (m && b && Q) { g->queueLength = g->maxMemory / g->blockSize;   }
  else if (m && b && q) { /* overspecified */ }

  // Isn't actually ever used anywhere, except for error checking below.
  g->maxMemory = g->queueLength * g->blockSize;

  if (g->queueLength < 128)
    sprintf(e, "ERROR: -queuelength must be at least 128.\n");
  if (g->blockSize == 0)
    sprintf(e, "ERROR: -blocksize must be more than zero.\n");
  if (g->maxMemory == 0)
    sprintf(e, "ERROR: -maxmemory must be more than zero.\n");

  if (g->inFiles.size() == 0)
    sprintf(e, "ERROR: no input-files supplied.\n");

  if (e.size() > 0) {
    fprintf(stderr, "usage: %s <input-files ...>\n", argv[0]);
    fprintf(stderr, "\n");
    for (char const *s : e)
      fputs(s, stderr);

    return(1);
  }

  ss->setLoaderBatchSize(8);
  ss->setLoaderQueueSize(1024);               //  There is no actual work done here;

  ss->setNumberOfWorkers(1);                  //  and we just dump all the data into

  ss->setWriterQueueSize(g->queueLength-64);  //  the output queue.
  ss->setInOrderOutput(true);

  ss->run(g, true);

  delete ss;
  delete g;

  return(0);
}
