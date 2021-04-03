
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
#include "system.H"

#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

#if defined(__FreeBSD__)
#include <stdlib.h>
#include <malloc_np.h>
#endif

#if defined(JEMALLOC)
#include "jemalloc/jemalloc.h"
#endif

#if !defined(__CYGWIN__) && !defined(_WIN32)
#include <sys/sysctl.h>
#endif



double
getTime(void) {
  struct timeval  tp;
  gettimeofday(&tp, NULL);
  return(tp.tv_sec + (double)tp.tv_usec / 1000000.0);
}



static
bool
getrusage(struct rusage &ru) {

  errno = 0;

  if (getrusage(RUSAGE_SELF, &ru) == -1) {
    fprintf(stderr, "getrusage(RUSAGE_SELF, ...) failed: %s\n",
            strerror(errno));
    return(false);
  }

  return(true);
}



static
bool
getrlimit(struct rlimit &rl) {

  errno = 0;

  if (getrlimit(RLIMIT_DATA, &rl) == -1) {
    fprintf(stderr, "getrlimit(RLIMIT_DATA, ...) failed: %s\n",
            strerror(errno));
    return(false);
  }

  return(true);
}



double
getCPUTime(void) {
  struct rusage  ru;
  double         tm = 0;

  if (getrusage(ru) == true)
    tm  = ((ru.ru_utime.tv_sec + ru.ru_utime.tv_usec / 1000000.0) +
           (ru.ru_stime.tv_sec + ru.ru_stime.tv_usec / 1000000.0));

  return(tm);
}



double
getProcessTime(void) {
  struct timeval tp;
  static double  st = 0.0;
  double         tm = 0;

  if (gettimeofday(&tp, NULL) == 0)
    tm  = tp.tv_sec + tp.tv_usec / 100000.0;

  if (st == 0.0)
    st = tm;

  return(tm - st);
}



uint64
getProcessSize(void) {
  struct rusage  ru;
  uint64         sz = 0;

  if (getrusage(ru) == true) {
    sz  = ru.ru_maxrss;
#ifndef __APPLE__     //  Everybody but MacOS returns kilobytes.
    sz *= 1024;       //  MacOS returns bytes.
#endif
  }

  return(sz);
}



uint64
getProcessSizeLimit(void) {
  struct rlimit rl;
  uint64        sz = ~uint64ZERO;

  if (getrlimit(rl) == true)
    sz = rl.rlim_cur;

  return(sz);
}



uint64
getBytesAllocated(void) {
  uint64 epoch     = 1;
  size_t epochLen  = sizeof(uint64);
  size_t active    = 0;
  size_t activeLen = sizeof(size_t);

#if defined(__FreeBSD__) || defined(JEMALLOC)

  mallctl("epoch", NULL, NULL, &epoch, epochLen);
  mallctl("stats.active", &active, &activeLen, NULL, 0);

#else

  active = getProcessSize();

#endif

  return(active);
}



#ifdef HW_PHYSMEM

//  MacOS, FreeBSD

uint64
getPhysicalMemorySize(void) {
  uint64  physMemory = 0;

  int     mib[2] = { CTL_HW, HW_PHYSMEM };
  size_t  len    = sizeof(uint64);

  errno = 0;

  if (sysctl(mib, 2, &physMemory, &len, NULL, 0) != 0)
    fprintf(stderr, "getPhysicalMemorySize()-- sysctl() failed to return CTL_HW, HW_PHYSMEM: %s\n", strerror(errno)), exit(1);

  if (len != sizeof(uint64)) {
#ifdef HW_MEMSIZE
    mib[1] = HW_MEMSIZE;
    len = sizeof(uint64);
    if (sysctl(mib, 2, &physMemory, &len, NULL, 0) != 0 || len != sizeof(uint64))
#endif
      fprintf(stderr, "getPhysicalMemorySize()-- sysctl() failed to return CTL_HW, HW_PHYSMEM: %s\n", strerror(errno)), exit(1);
  }

  return(physMemory);
}

#else

//  Linux, FreeBSD

uint64
getPhysicalMemorySize(void) {
  uint64  physPages  = sysconf(_SC_PHYS_PAGES);
  uint64  pageSize   = sysconf(_SC_PAGESIZE);
  uint64  physMemory = physPages * pageSize;

  return(physMemory);
}

#endif




//  Return the size of a page of memory.  Every OS we care about (MacOS, FreeBSD, Linux)
//  claims to have getpagesize().
//
uint64
getPageSize(void) {
  return(getpagesize());
}
