
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


void  AS_UTL_catchCrash(int sig_num, siginfo_t *info, void *ctx);

void  AS_UTL_installCrashCatcher(const char *filename);



#endif  //  SYSTEM_H
