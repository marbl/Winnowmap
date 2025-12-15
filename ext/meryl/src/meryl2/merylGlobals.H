
/******************************************************************************
 *
 *  This file is part of meryl, a genomic k-kmer counter with nice features.
 *
 *  This software is based on:
 *    'Canu' v2.0              (https://github.com/marbl/canu)
 *  which is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "types.H"
#include "system.H"

//
//  A global class to hold:
//   - maximum memory usage allowed
//   - maximum CPU usage allowed
//   - verbosity level
//
//  It also handles some debug-level command line functions that don't really
//  fit anywhere else.
//

class merylGlobals {
public:
  merylGlobals()  {};
  ~merylGlobals() {};

  int    initialize(int argc, char **argv) {
    int a = AS_configure(argc, argv);

    _aM = getMaxMemoryAllowed();   //  Remember any system-imposed
    _aT = getMaxThreadsAllowed();  //  limit on cpu/threads.

    return a;
  };

  bool    processDebugOption (int &arg, char **argv, std::vector<char const *> &err);
  bool    processGlobalOption(int &arg, char **argv, std::vector<char const *> &err);

  void    beQuiet(void)               { _v = 0; }
  void    increaseVerbosity(void)     { _v++;   }

  bool    showNothing     (void)      { return(_v == 0); }
  bool    showStandard    (void)      { return(_v >= 1); }
  bool    showProgress    (void)      { return(_v >= 2); }
  bool    showConstruction(void)      { return(_v >= 3); }
  bool    showDetails     (void)      { return(_v >= 4); }
  bool    showEverything  (void)      { return(_v >= 5); }

  void    enableProgressReport(void)  { _p = true; }
  bool    showProgressReport(void)    { return _p; }

  void    enableConfigureOnly(void)   { _c = true; }
  bool    stopAfterConfigure(void)    { return _c; }

  uint64 &allowedMemory(void)         { return _aM; }
  uint32 &allowedThreads(void)        { return _aT; }

private:
  uint32 _v  = 1;       //  Verbosity level, 'standard'.
  bool   _p  = false;   //  Show progress, does nothing yet.

  //nt32 _k  = 0;       //  Both k-mer size and label size are
  //nt32 _l  = 0;       //  stored in the kmer class.

  bool   _c  = false;   //  Stop after configure.

  uint64 _aM = 0;       //  Max amount of memory to use.
  uint32 _aT = 0;       //  Max number of threads to use.
};

extern merylGlobals globals;

