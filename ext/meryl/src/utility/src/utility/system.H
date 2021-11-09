
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

#ifndef SYSTEM_H
#define SYSTEM_H

#include "types.H"

#include <signal.h>



double   getTime(void);

double   getCPUTime(void);
double   getProcessTime(void);

uint64   getProcessSize(void);
uint64   getProcessSizeLimit(void);

uint64   getBytesAllocated(void);

uint64   getPhysicalMemorySize(void);

uint64   getPageSize(void);

uint64   getMaxMemoryAllowed(void);

//
//  getMaxThreadsAllowed() returns the minimum value of 'limit' and any of a
//  variety of environment variables specifying a thread limit.  It does NOT
//  configure OpenMP to use that many threads.
//
//  AS_configure() (runtime.H) calls:
//    setNumThreads(getMaxThreadsAllowed(limit))
//  User code is then free to parse args and call setNumThreads() to change
//  the number of threads to whatever is specified on the command line.
//
//  getNumThreads() and getNumThreadsActive() return the current limit on the
//  number of threads and the current number of threads running, respectively.
//
//  setNumThreads() will configure OpenMP to use at most the supplied number
//  of threads.  It returns this number.
//
//  getThreadNum() returns 0 .. numThreadsActive-1 for each thread.
//
uint32   getMaxThreadsAllowed(uint32 limit=uint32max);
uint32   getNumThreads(void);
uint32   getNumThreadsActive(void);

uint32   getThreadNum(void);

uint32   setNumThreads(char const *opt);
uint32   setNumThreads(uint32 thr);

void  AS_UTL_catchCrash(int sig_num, siginfo_t *info, void *ctx);

void  AS_UTL_installCrashCatcher(const char *filename);



#endif  //  SYSTEM_H
