
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
    tm  = tp.tv_sec + tp.tv_usec / 1000000.0;

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
  uint64        sz = uint64max;

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



uint64
getPhysicalMemorySize(void) {
  uint64  physPages  = sysconf(_SC_PHYS_PAGES);
  uint64  pageSize   = sysconf(_SC_PAGESIZE);
  uint64  physMemory = physPages * pageSize;

  return(physMemory);
}



//  Return the size of a page of memory.  Every OS we care about (MacOS,
//  FreeBSD, Linux) claims to have getpagesize().
//
uint64
getPageSize(void) {
  return(getpagesize());
}



//  Query the machine or the environment to find any memory size limit.  If
//  there is no environment limit, the physical memory size is returned.
//
//  Slurm variables (from sbatch man page).
//    SLURM_MEM_PER_CPU
//      Set if --mem-per-cpu is supplied to sbatch.
//      "SLURM_MEM_PER_CPU=2048" for a request of --mem-per-cpu=2g
//
//    SLURM_MEM_PER_NODE
//      Set if --mem is supplied to sbatch.
//      "SLURM_MEM_PER_NODE=5120" for a request of --mem=5g
//
//    SLURM_MEM_PER_GPU
//      Requested memory per allocated GPU.
//        Only set if the --mem-per-gpu option is specified.
//        Not checked for below.
//
//  There doesn't appear to be a comparable environment variable for SGE.
//
//  PBS/OpenPBS/PBS Pro variables.
//    PBS_RESC_MEM
//    TORQUE_RESC_MEM  (probably obsolete)
//      Potentially memory in bytes.
//
//
uint64
getMaxMemoryAllowed(void) {
  char    *env, *cpu;
  uint64   maxmem = getPhysicalMemorySize();

  cpu = getenv("SLURM_JOB_CPUS_PER_NODE");
  env = getenv("SLURM_MEM_PER_CPU");
  if (env && cpu)
    maxmem = strtouint64(cpu) * strtouint64(env) * 1024 * 1024;

  env = getenv("SLURM_MEM_PER_NODE");
  if (env)
    maxmem = strtouint64(env) * 1024 * 1024;

  env = getenv("PBS_RESC_MEM");
  if (env)
    maxmem = strtouint64(env);

  return(maxmem);
}



//  There is a bit of a race condition in here.  On our grid, at least, a
//  multi-cpu interactive job sets both SLURM_JOB_CPUS_PER_NODE and
//  OMP_NUM_THREADS - but sets the former to the correct value and the
//  latter to one.
//
//  Because of this, we let the grid variables overwrite the OpenMP
//  variable, and further reset OpenMP to use whatever the grid has
//  told us to use.
//
//  OpenMP variables.
//    OMP_NUM_THREADS
//     - we don't query this, and instead use omp_get_max_threads(),
//       because if OMP_NUM_THREADS isn't set, the function will
//       return the number of CPUs on the host.
//
//  Slurm variables (from sbatch man page).
//    SLURM_CPUS_ON_NODE
//     - Number of CPUS on the allocated node.
//
//    SLURM_JOB_CPUS_PER_NODE
//     - --cpus-per-node
//     - Count of processors available to the job on this node. Note the
//       select/linear plugin allocates entire nodes to jobs, so the value
//       indicates the total count of CPUs on the node. The select/cons_res
//       plugin allocates individual processors to jobs, so this number
//       indicates the number of processors on this node allocated to the
//       job.
//
//    SLURM_JOB_NUM_NODES
//     - total number of nodes in the job's resource allocation
//
//  PBS/OpenPBS/PBS Pro variables (from Torque 9.0.3).
//    PBS_NUM_NODES    - Number of nodes allocated to the job
//    PBS_NUM_PPN      - Number of procs per node allocated to the job
//    PBS_NCPUS        - (older version of PBS_NUM_PPN?)
//    PBS_NP           - Number of execution slots (cores) for the job
//    TORQUE_RESC_PROC - (can't find any doc on this)
//
//  SGE variables.
//    NSLOTS
//
uint32
getMaxThreadsAllowed(uint32 limit) {
  char    *env;
  uint32   nAllowed = omp_get_max_threads();

  env = getenv("SLURM_JOB_CPUS_PER_NODE");
  if (env)
    nAllowed = strtouint32(env);

  env = getenv("PBS_NCPUS");
  if (env)
    nAllowed = strtouint32(env);

  env = getenv("PBS_NUM_PPN");
  if (env)
    nAllowed = strtouint32(env);

  env = getenv("NSLOTS");
  if (env)
    nAllowed = strtouint32(env);

  nAllowed = std::min(limit, nAllowed);

  return(nAllowed);
}



uint32
getNumThreads(void) {
  return(omp_get_max_threads());
}



uint32
getNumThreadsActive(void) {
  return(omp_get_num_threads());
}



uint32
getThreadNum(void) {
  return(omp_get_thread_num());
}



uint32
setNumThreads(char const *opt) {
  return(setNumThreads(strtouint32(opt)));
}



uint32
setNumThreads(uint32 thr) {
  omp_set_num_threads(thr);
  return(omp_get_max_threads());
}
