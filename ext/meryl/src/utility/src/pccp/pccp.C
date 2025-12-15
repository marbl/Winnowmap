#include "types.H"
#include "files.H"
#include "system.H"
#include "math.H"

#include <vector>


class cpBuf {
public:
  cpBuf(FILE *F, uint32 max) {
    _b    = new uint8 [max];
    _bLen = merylutil::loadFromFile(_b, "stuff", max, F, false);
  }

  ~cpBuf() {
    delete [] _b;
  }

  uint32  _bLen = 0;
  uint8  *_b    = nullptr;
};



class cpBufState {
public:
  cpBufState(char const *inname, char const *otpath) {
    snprintf(inPath, FILENAME_MAX, "%s", inname);
    snprintf(otPath, FILENAME_MAX, "%s/%s", otpath, basename(inname));

    fileSize = merylutil::sizeOfFile(inPath);

    if ((merylutil::fileExists(otPath) == true) &&
        (merylutil::sizeOfFile(otPath) != fileSize))
      merylutil::unlink(otPath);

    if (merylutil::fileExists(otPath) == false) {
      inFile = merylutil::openInputFile(inPath);
      otFile = merylutil::openOutputFile(otPath);
    }
  };

  //  If the output was opened, close the files and update the timestamp.
  //  If it wasn't opened, it existed already and we didn't do the copy.
  ~cpBufState() {
    if (otFile) {
      merylutil::muFileTime times;

      merylutil::closeFile(inFile);
      merylutil::closeFile(otFile);

      times.getTimeOfFile(inPath);
      times.setTimeOfFile(otPath);
    }
  };

  char const *basename(char const *name) {
    uint32 l = strlen(name);

    for (uint32 ii=l; --ii; ) {
      if (name[ii] == '/')
        return(name + ii + 1);
    }

    return(name);
  }

  cpBuf      *read(void) {
    cpBuf  *b = new cpBuf(inFile, 1048576);
    double  n = getTime();

    inStart = std::min(inStart, n);
    inEnd   = std::max(inEnd,   n);

    if (b->_bLen == 0) {
      inEnd = getTime();
      delete b;
      return(nullptr);
    }

    return(b);
  };

  void        write(cpBuf *b) {
    double  n = getTime();

    otStart = std::min(otStart, n);
    otEnd   = std::max(otEnd,   n);

    merylutil::writeToFile(b->_b, "stuff", b->_bLen, otFile);

    delete b;
  }

  char        inPath[FILENAME_MAX];
  char        otPath[FILENAME_MAX];

  FILE       *inFile = nullptr;
  FILE       *otFile = nullptr;

  uint64      bufMax   = 1048576;
  uint64      fileSize = 0;

  double      inStart=DBL_MAX, inEnd=0;
  double      otStart=DBL_MAX, otEnd=0;

  merylutil::md5sum  md5;
};




void *bufReader(void *G) {
  cpBufState *g = (cpBufState *)G;
  return(g->read());
}

void  bufWorker(void *G, void *T, void *S) {
  cpBufState *g = (cpBufState *)G;
  cpBuf      *s = (cpBuf      *)S;

  g->md5.addBlock(s->_b, s->_bLen);
}

void  bufWriter(void *G, void *S) {
  cpBufState *g = (cpBufState *)G;
  cpBuf      *s = (cpBuf      *)S;
  g->write(s);
}

void  bufStatus(void *G, uint64 numberLoaded, uint64 numberComputed, uint64 numberOutput) {
  cpBufState *g = (cpBufState *)G;

  double  thisTime = getTime();

  double  inSize  = numberLoaded * g->bufMax;
  double  inPerc  = 100.0 * inSize / g->fileSize;
  double  inSpeed =         inSize / (g->inEnd - g->inStart);

  double  otSize  = numberOutput * g->bufMax;
  double  otPerc  = 100.0 * otSize / g->fileSize;
  double  otSpeed =         otSize / (g->otEnd - g->otStart);

  if (inPerc > 100)   inPerc = 100.0;
  if (otPerc > 100)   otPerc = 100.0;

  if (g->otFile) {
    fprintf(stderr, "\033[3A");
    fprintf(stderr, "  INPUT:  %6.2f%%  %8.2f MB  %6.2f MB/sec\n", inPerc, inSize / 1048576.0, inSpeed / 1048576.0);
    fprintf(stderr, "  BUFFER:          %8.2f MB\n", (inSize - otSize) / 1048576.0);
    fprintf(stderr, "  OUTPUT: %6.2f%%  %8.2f MB  %6.2f MB/sec\n", otPerc, otSize / 1048576.0, otSpeed / 1048576.0);
  }
}


int
main(int argc, char **argv) {
  std::vector<char const *>   errors;
  std::vector<char const *>   infiles;
  char const                 *otpath  = nullptr;

  uint32  lqSize = 128;    //  Loader Queue size
  uint32  wqSize = 16384;  //  Writer Queue size

  for (int arg=1; arg<argc; arg++) {
    if      (strcmp(argv[arg], "-b") == 0)
      wqSize = strtouint32(argv[++arg]);

    else if (merylutil::fileExists(argv[arg]))
      infiles.push_back(argv[arg]);

    else if (merylutil::directoryExists(argv[arg]) && (arg == argc-1))
      otpath = argv[arg];

    else if (merylutil::directoryExists(argv[arg]))
      ;

    else
      sprintf(errors, "ERROR: '%s' is neither an input file nor an output directory.\n", argv[arg]);
  }

  if (infiles.size() == 0)
    sprintf(errors, "ERROR: no input-files supplied.\n");
  if (otpath == nullptr)
    sprintf(errors, "ERROR: no output-directory supplied.\n");

  if (errors.size() > 0) {
    fprintf(stderr, "usage: %s [-b s] <input-files ...> <output-directory>\n", argv[0]);
    fprintf(stderr, "  -b s     limit to 's' 1-MB buffers\n");
    for (char const *e : errors)
      fputs(e, stderr);

    return(1);
  }

  for (char const *infile : infiles) {
    cpBufState *g  = new cpBufState(infile, otpath);
    sweatShop  *ss = new sweatShop(bufReader, bufWorker, bufWriter, bufStatus);

    ss->setLoaderQueueSize(lqSize);   //  Just comuting the md5, so we don't need lots of input buffering.
    ss->setNumberOfWorkers(1);        //
    ss->setWriterQueueSize(wqSize);   //  But if we wait for writing, keep filling the queue.
    ss->setInOrderOutput(true);

    fprintf(stderr, "%s -> %s\n", infile, otpath);

    if (g->otFile) {
      fprintf(stderr, "\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "\n");

      ss->run(g, true);

      g->md5.finalize();

      fprintf(stderr, "\033[3A");
      fprintf(stderr, "  MD5:    %s\n", g->md5.toString());
      fprintf(stderr, "\033[3B");
    }

    delete ss;
    delete g;
  }

  return(0);
}
