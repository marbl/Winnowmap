
/******************************************************************************
 *
 *  This file is part of meryl-utility, a collection of miscellaneous code
 *  used by Meryl, Canu and others.
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
#include "files.H"
#include "system.H"

#include "version.H"

#ifdef X86_GCC_LINUX
#include <fpu_control.h>
#endif

#ifdef _GLIBCXX_PARALLEL
#include <parallel/algorithm>
#include <parallel/settings.h>
#endif



//  Set the x86 FPU control word to force double precision rounding
//  rather than `extended' precision rounding. This causes base
//  calls and quality values on x86 GCC-Linux (tested on RedHat
//  Linux) machines to be identical to those on IEEE conforming UNIX
//  machines.

#ifdef X86_GCC_LINUX

void
setFPU(void) {
  fpu_control_t fpu_cw = ( _FPU_DEFAULT & ~_FPU_EXTENDED ) | _FPU_DOUBLE;

  _FPU_SETCW( fpu_cw );
}

#else

void
setFPU(void) {
}

#endif



//  The sizes that GNU decides to enable parallelization for are
//  ludicrously small.  A benchmark (Dec 2019, gcc6) showed that a
//  sort-heavy code was a whopping 10% faster (wall clock 61 vs 69 seconds)
//  but used 10x the CPU time (707 vs 69 seconds).
//
//  Further, the sort is no longer in-place.  This matters, significantly,
//  for our overlap sorting, since we fill memory with as many overlaps as
//  possible, and if we don't sort in-place, we run out of memory.
//
//  So, if the silly user has enabled parallelization, turn it off for
//  sorting by setting the minimal size to something large.

#ifdef _GLIBCXX_PARALLEL_SETTINGS_H

void
setSequentialSorting(void) {
  __gnu_parallel::_Settings s = __gnu_parallel::_Settings::get();

  s.sort_minimal_n = UINT64_MAX;

  __gnu_parallel::_Settings::set(s);
}

#else

void
setSequentialSorting(void) {
}

#endif





int
AS_configure(int argc, char **argv, uint32 maxThreads) {

  setFPU();
  setSequentialSorting();

  AS_UTL_installCrashCatcher(argv[0]);

  getProcessTime();   //  To set the process start time.

  //  Enable threads to whatever limit is given to us or whatever the
  //  envrionment says.

  setNumThreads(getMaxThreadsAllowed(maxThreads));

  //  Parse any common options from the command line.  AS_configure() was
  //  designed to consume options, but this isn't needed at the moment.

  for (int32 i=0; i<argc; i++) {
    if (strcmp(argv[i], "--version") == 0) {
      fputs(MERYL_UTILITY_VERSION, stderr);
      exit(0);
    }
  }

  return(argc);
}
