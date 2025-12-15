
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

#ifndef SWEATSHOP_H
#define SWEATSHOP_H

#include "types.H"

#include <pthread.h>

class sweatShopWorker;
class sweatShopState;

class sweatShop {
public:
  sweatShop(void*(*loaderfcn)(void *G),
            void (*workerfcn)(void *G, void *T, void *S),
            void (*writerfcn)(void *G, void *S),
            void (*statusfcn)(void *G, uint64 numberLoaded, uint64 numberComputed, uint64 numberOutput) = nullptr);
  ~sweatShop();

  void        setThreadData(uint32 t, void *x);

  void        setNumberOfWorkers(uint32 x)         { _numberOfWorkers = x; }

  void        setLoaderBatchSize(uint32 batchSize) { _loaderBatchSize = batchSize; }
  void        setLoaderQueueSize(uint32 queueSize) { _loaderQueueSize = queueSize;  _loaderQueueMax = queueSize; }

  void        setWorkerBatchSize(uint32 batchSize) { _workerBatchSize = batchSize; }

  void        setWriterQueueSize(uint32 queueSize) { _writerQueueSize = queueSize;  _writerQueueMax = queueSize; }

  void        setInOrderOutput(bool o)             { _writeInOrder = o; }

  void        run(void *user=0L, bool beVerbose=false);
private:

  //  Stubs that forward control from the c-based pthread to the
  //  corresponding methods in this class, then those methods forward the
  //  call to the user supplied functions that do the actual work.

  friend void *_sload(void *ss);   void *loader(void);
  friend void *_swork(void *ss);   void *worker(sweatShopWorker *workerData);
  friend void *_swrit(void *ss);   void *writer(void);
  friend void *_sstat(void *ss);   void *status(void);

  void *(*_userLoader)(void *global);
  void  (*_userWorker)(void *global, void *thread, void *thing);
  void  (*_userWriter)(void *global, void *thing);
  void  (*_userStatus)(void *global, uint64 numberLoaded, uint64 numberComputed, uint64 numberOutput);


  //  Utilities for the loader thread
  void    loaderAddToLocal(sweatShopState *&tail, sweatShopState *&head, sweatShopState *thisState);
  void    loaderAppendToGlobal(sweatShopState *&tail, sweatShopState *&head, uint32 num);

  //  Utilities for the writer thread
  void    writerWrite(sweatShopState *w);

  pthread_mutex_t        _stateMutex;


  void                  *_globalUserData;

  sweatShopState        *_writerP;  //  Where output takes stuff from, the tail
  sweatShopState        *_workerP;  //  Where computes happen, the middle
  sweatShopState        *_loaderP;  //  Where input is put, the head

  bool                   _showStatus;
  bool                   _writeInOrder;

  uint32                 _loaderQueueSize, _loaderQueueMax;
  uint32                 _loaderBatchSize;
  uint32                 _workerBatchSize;
  uint32                 _writerQueueSize, _writerQueueMax;

  uint32                 _numberOfWorkers;

  sweatShopWorker       *_workerData;

  uint64                 _numberLoaded;
  uint64                 _numberComputed;
  uint64                 _numberOutput;
};

#endif  //  SWEATSHOP_H
