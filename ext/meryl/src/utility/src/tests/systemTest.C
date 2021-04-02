
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

#include "system.H"


int
main(int argc, char **argv) {

  fprintf(stderr, "This is most useful if you run it as 'time %s'\n", argv[0]);
  fprintf(stderr, "then compare the two reports.\n");
  fprintf(stderr, "\n");

  getProcessTime();  //  initialize.

  double   *array = new double [1024 * 1024];

  for (int jj=0; jj<128; jj++)
    for (int ii=0; ii<1024*1024; ii++)
      array[ii] = sin(ii) + cos(jj);

  double     t = getTime();
  time_t     T = t;

  fprintf(stderr, "getTime()               %f\n",     t);
  fprintf(stderr, "getTime()  (ctime)      %s\n",     ctime(&T));
  fprintf(stderr, "\n");
  fprintf(stderr, "getCPUTime              %f\n",     getCPUTime());
  fprintf(stderr, "getProcessTime()        %f\n",     getProcessTime());
  fprintf(stderr, "getProcessSize()        %lu\n",    getProcessSize());
  fprintf(stderr, "getProcessSizeLimit()   %lu\n",    getProcessSizeLimit());
  fprintf(stderr, "getProcessSizeLimit()   %lu GB\n", getProcessSizeLimit() >> 30);
  fprintf(stderr, "getBytesAllocated()     %lu\n",    getBytesAllocated());
  fprintf(stderr, "\n");
  fprintf(stderr, "getPhysicalMemorySize() %lu\n",    getPhysicalMemorySize());
  fprintf(stderr, "getPhysicalMemorySize() %lu GB\n", getPhysicalMemorySize() >> 30);
  fprintf(stderr, "getPageSize()           %lu\n",    getPageSize());
  fprintf(stderr, "\n");
  fprintf(stderr, "getMaxMemoryAllowed()   %lu\n",    getMaxMemoryAllowed());
  fprintf(stderr, "getMaxMemoryAllowed()   %lu GB\n", getMaxMemoryAllowed() >> 30);
  fprintf(stderr, "getMaxThreadsAllowed()  %u\n",     getMaxThreadsAllowed());
  fprintf(stderr, "\n");

  return(0);
}
