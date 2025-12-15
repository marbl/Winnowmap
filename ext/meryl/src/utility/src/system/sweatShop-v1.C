
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
#include "sweatShop-v1.H"

#include <pthread.h>
#include <sched.h>  //  pthread scheduling stuff

class sweatShopWorker {                          //  Description of a single worker
public:                                          //  thread/queue.
  pthread_t         threadID;

  uint32            numComputed    = 0;

  sweatShop        *sweatshop      = nullptr;
  void             *threadUserData = nullptr;

  sweatShopState  **workerQueue    = nullptr;
  uint32            workerQueueLen = 0;
};


class sweatShopState {                           //  Object holding a single element for
public:                                          //  load, compute and write.
  sweatShopState(void *userData) : _userdata(userData) {};

  void             *_userdata  = nullptr;        //  User supplied data element.
  bool              _computed  = false;          //  Is this object computed?
  bool              _outputted = false;          //  Is this object written?
  sweatShopState   *_next      = nullptr;        //  Pointer to any next object.
};


//  Stubs that forward control from pthreads back to the sweatShop.
void *_sload(void *ss_) { return ((sweatShop       *)ss_)->loader(); }
void *_swork(void *sw_) { return ((sweatShopWorker *)sw_)->sweatshop->worker((sweatShopWorker *)sw_); }
void *_swrit(void *ss_) { return ((sweatShop       *)ss_)->writer(); }
void *_sstat(void *ss_) { return ((sweatShop       *)ss_)->status(); }


sweatShop::sweatShop(void*(*loaderfcn)(void *G),
                     void (*workerfcn)(void *G, void *T, void *S),
                     void (*writerfcn)(void *G, void *S),
                     void (*statusfcn)(void *G, uint64 numberLoaded, uint64 numberComputed, uint64 numberOutput)) {

  _userLoader       = loaderfcn;
  _userWorker       = workerfcn;
  _userWriter       = writerfcn;
  _userStatus       = statusfcn;

  _globalUserData   = nullptr;

  _writerP          = nullptr;
  _workerP          = nullptr;
  _loaderP          = nullptr;

  _showStatus       = false;
  _writeInOrder     = true;

  _loaderQueueSize  = 1024;
  _loaderQueueMax   = 10240;
  _loaderBatchSize  = 1;
  _workerBatchSize  = 1;
  _writerQueueSize  = 4096;

  _numberOfWorkers  = 2;

  _workerData       = nullptr;

  _numberLoaded     = 0;
  _numberComputed   = 0;
  _numberOutput     = 0;
}


sweatShop::~sweatShop() {
  delete [] _workerData;
}



void
sweatShop::setThreadData(uint32 t, void *x) {
  if (_workerData == nullptr)
    _workerData = new sweatShopWorker [_numberOfWorkers];

  if (t >= _numberOfWorkers)
    fprintf(stderr, "sweatShop::setThreadData()-- worker ID " F_U32 " more than number of workers=" F_U32 "\n", t, _numberOfWorkers), exit(1);

  _workerData[t].threadUserData = x;
}



//  Build a list of states to add in one swoop
//
void
sweatShop::loaderAddToLocal(sweatShopState *&tail, sweatShopState *&head, sweatShopState *thisState) {

  thisState->_next  = nullptr;

  if (tail) {
    head->_next = thisState;
    head        = thisState;
  } else {
    tail = head = thisState;
  }
}


//  Add a bunch of new states to the queue.
//
void
sweatShop::loaderAppendToGlobal(sweatShopState *&tail, sweatShopState *&head, uint32 num) {
  int err;

  if ((tail == nullptr) || (head == nullptr))
    return;

  err = pthread_mutex_lock(&_stateMutex);
  if (err != 0)
    fprintf(stderr, "sweatShop::loaderAppend()--  Failed to lock mutex (%d).  Fail.\n", err), exit(1);

  if (_loaderP == nullptr) {
    _writerP      = tail;
    _workerP      = tail;
    _loaderP      = head;
  } else {
    _loaderP->_next = tail;
  }
  _loaderP        = head;

  _numberLoaded += num;

  err = pthread_mutex_unlock(&_stateMutex);
  if (err != 0)
    fprintf(stderr, "sweatShop::loaderAppend()--  Failed to unlock mutex (%d).  Fail.\n", err), exit(1);

  tail = nullptr;
  head = nullptr;
}



void*
sweatShop::loader(void) {
  struct timespec   naptime;
  sweatShopState   *tail       = nullptr;  //  A local list, to reduce the number of times we
  sweatShopState   *head       = nullptr;  //  lock the global list.
  uint32            numLoaded  = 0;

  naptime.tv_sec      = 0;
  naptime.tv_nsec     = 166666666ULL;  //  1/6 second

  while (1) {
    void *object = NULL;

    while (_numberLoaded > _numberComputed + _loaderQueueSize)  //  Sleep if the queue is too big.
      nanosleep(&naptime, nullptr);

    //  If a userLoader function exists, use it to load the data object, then
    //  make a new state for that object.

    if (_userLoader)
      object = (*_userLoader)(_globalUserData);

    sweatShopState  *thisState = new sweatShopState(object);

    //  If there is no user pointer, we've run out of inputs.
    //  Push on the empty state to the local list, force an append
    //  to the global list, and exit this loader function.

    if (thisState->_userdata == nullptr) {
      loaderAddToLocal(tail, head, thisState);
      loaderAppendToGlobal(tail, head, numLoaded + 1);

      return(nullptr);
    }

    //  Otherwise, we've loaded a user object.  Push it onto the local list,
    //  then merge into the global list if the local list is long enough.

    loaderAddToLocal(tail, head, thisState);
    numLoaded++;

    if (numLoaded >= _loaderBatchSize) {
      loaderAppendToGlobal(tail, head, numLoaded);
      numLoaded = 0;
    }
  }

  return(nullptr);  //  Never returns.
}



void*
sweatShop::worker(sweatShopWorker *workerData) {

  struct timespec   naptime;
  naptime.tv_sec      = 0;
  naptime.tv_nsec     = 50000000ULL;

  bool    moreToCompute = true;
  int     err;

  while (moreToCompute) {

    //  Usually beacuse some worker is taking a long time, and the
    //  output queue isn't big enough.
    //
    while (_numberOutput + _writerQueueSize < _numberComputed)
      nanosleep(&naptime, nullptr);

    //  Grab the next state.  We don't grab it if it's the last in the
    //  queue (else we would fall off the end) UNLESS it really is the
    //  last one.
    //
    err = pthread_mutex_lock(&_stateMutex);
    if (err != 0)
      fprintf(stderr, "sweatShop::worker()--  Failed to lock mutex (%d).  Fail.\n", err), exit(1);

    for (workerData->workerQueueLen = 0; ((workerData->workerQueueLen < _workerBatchSize) &&
                                          (_workerP) &&
                                          ((_workerP->_next != nullptr) || (_workerP->_userdata == nullptr))); workerData->workerQueueLen++) {
      workerData->workerQueue[workerData->workerQueueLen] = _workerP;
      _workerP = _workerP->_next;
    }

    if (_workerP == nullptr)
      moreToCompute = false;

    err = pthread_mutex_unlock(&_stateMutex);
    if (err != 0)
      fprintf(stderr, "sweatShop::worker()--  Failed to lock mutex (%d).  Fail.\n", err), exit(1);


    if (workerData->workerQueueLen == 0) {
      //  No work, sleep a bit to prevent thrashing the mutex and resume.
      nanosleep(&naptime, nullptr);
      continue;
    }

    //  Execute
    //
    for (uint32 x=0; x<workerData->workerQueueLen; x++) {
      sweatShopState *ts = workerData->workerQueue[x];

      if (ts && ts->_userdata) {
        if (_userWorker)
          (*_userWorker)(_globalUserData, workerData->threadUserData, ts->_userdata);
        ts->_computed = true;
        workerData->numComputed++;
      } else {
        //  When we really do run out of stuff to do, we'll end up here
        //  (only one thread will end up in the other case, with
        //  something to do and moreToCompute=false).  If it's actually
        //  the end, skip the sleep and just get outta here.
        //
        if (moreToCompute == true) {
          fprintf(stderr, "WARNING!  Worker is sleeping because the reader is slow!\n");
          nanosleep(&naptime, nullptr);
        }
      }
    }
  }

  //fprintf(stderr, "sweatShop::worker exits.\n");
  return(nullptr);
}


void
sweatShop::writerWrite(sweatShopState *w) {

  if (_userWriter)
    (*_userWriter)(_globalUserData, w->_userdata);
  _numberOutput++;

  w->_outputted = true;
}


void*
sweatShop::writer(void) {
  sweatShopState  *deleteState = nullptr;
  struct timespec naptime1 = { .tv_sec = 0, .tv_nsec = 5000000ULL };
  struct timespec naptime2 = { .tv_sec = 0, .tv_nsec = 5000000ULL };


  while ((_writerP        != nullptr) &&
         (_writerP->_userdata != nullptr)) {

    //  If a complete result, write it.
    if ((_writerP->_computed  == true) &&
        (_writerP->_outputted == false)) {
      writerWrite(_writerP);
      continue;
    }

    //  If we can write output out-of-order, search ahead
    //  for any results and output them.
    //  if (_outOfOrder == true)
    if (_writeInOrder == false) {
      for (sweatShopState *ss = _writerP; ss != nullptr; ss = ss->_next)
        if ((ss->_computed  == true) &&
            (ss->_outputted == false)) {
          writerWrite(ss);
        }
    }

    //  If no next, wait for input to appear.  We can't purge this node
    //  from the list until there is a next, else we lose the list!
    if (_writerP->_next == nullptr) {
      nanosleep(&naptime1, nullptr);
      continue;
    }

    //  If already output, remove the node.
    if (_writerP->_outputted == true) {
      sweatShopState *ds = _writerP;
      _writerP           = _writerP->_next;

      delete ds;
      continue;
    }

    //  Otherwise, we need to wait for a state to appear on the queue.
    nanosleep(&naptime2, nullptr);
  }

  //  Tell status to stop.
  _writerP = nullptr;

  return(nullptr);
}


//  This thread not only shows a status message, but it also updates the critical shared variable
//  _numberComputed.  Worker threads use this to throttle themselves.  Thus, even if _showStatus is
//  not set, and this thread doesn't _appear_ to be doing anything useful....it is.
//

void
sweatShopStatus(double startTime, uint64 numberLoaded, uint64 numberComputed, uint64 numberOutput) {
  double thisTime  = getTime();
  uint64 deltaOut  = 0;
  uint64 deltaCPU  = 0;
  double cpuPerSec = 0;

  if (numberComputed > numberOutput)
    deltaOut = numberComputed - numberOutput;
  if (numberLoaded > numberComputed)
    deltaCPU = numberLoaded - numberComputed;

  cpuPerSec = numberComputed / (thisTime - startTime);

  fprintf(stderr, " %6.1f/s - %8" F_U64P " loaded; %8" F_U64P " queued for compute; %8" F_U64P " finished; %8" F_U64P " written; %8" F_U64P " queued for output)\r",
          cpuPerSec, numberLoaded, deltaCPU, numberComputed, numberOutput, deltaOut);
}



void*
sweatShop::status(void) {
  struct timespec   naptime;
  naptime.tv_sec      = 0;
  naptime.tv_nsec     = 250000000ULL;

  double  startTime  = getTime() - 0.001;
  uint64  readjustAt = 16384;

  while (_writerP) {
    uint32 nc = 0;
    for (uint32 i=0; i<_numberOfWorkers; i++)
      nc += _workerData[i].numComputed;
    _numberComputed = nc;

    if (_showStatus) {
      if (_userStatus)
        (*_userStatus)(_globalUserData, _numberLoaded, _numberComputed, _numberOutput);
      else
        sweatShopStatus(startTime, _numberLoaded, _numberComputed, _numberOutput);

      fflush(stderr);
    }

    if (_numberComputed > readjustAt) {                               //  Adjust input queue size based on
      double cpuPerSec = _numberComputed / (getTime() - startTime);   //  current performance, but don't let
                                                                      //  it get too big (above user supplied
      readjustAt       += (uint64)(2 * cpuPerSec);                    //  limit) or too small (below number
      _loaderQueueSize  = (uint32)(5 * cpuPerSec);                    //  of compute threads).

      _loaderQueueSize = std::min(_loaderQueueSize, _loaderQueueMax);       //  Too big!
      _loaderQueueSize = std::max(_loaderQueueSize, _numberOfWorkers * 2);  //  Too small!
    }

    nanosleep(&naptime, nullptr);
  }

  if (_showStatus) {
    if (_userStatus)
      (*_userStatus)(_globalUserData, _numberLoaded, _numberComputed, _numberOutput);
    else
      sweatShopStatus(startTime, _numberLoaded, _numberComputed, _numberOutput);

    fprintf(stderr, "\n");
    fflush(stderr);
  }

  //fprintf(stderr, "sweatShop::status exits.\n");
  return(nullptr);
}



static
void
_rE(int32 err, char const *msg=nullptr, uint32 idx=uint32max) {
  if (err) {
    if (idx == uint32max)   fprintf(stderr, msg,      strerror(err));
    else                    fprintf(stderr, msg, idx, strerror(err));
    exit(1);
  }
}


void
sweatShop::run(void *user, bool beVerbose) {
  struct timespec     naptime = { 0, 250000 };

  pthread_attr_t      threadAttr;
  pthread_t           threadIDloader;
  pthread_t           threadIDwriter;
  pthread_t           threadIDstats;

  _globalUserData = user;
  _showStatus     = beVerbose;

  //  Configure everything ahead of time.

  if (_workerBatchSize < 1)
    _workerBatchSize = 1;

  if (_workerData == nullptr)
    _workerData = new sweatShopWorker [_numberOfWorkers];

  for (uint32 i=0; i<_numberOfWorkers; i++) {
    _workerData[i].sweatshop   = this;
    _workerData[i].workerQueue = new sweatShopState * [_workerBatchSize];
  }

  //  Open the doors.

  _rE(pthread_mutex_init(&_stateMutex, NULL),                            "sweatShop::run()--  Failed to configure pthreads (state mutex): %s.\n");
  _rE(pthread_attr_init(&threadAttr),                                    "sweatShop::run()--  Failed to configure pthreads (attr init): %s.\n");
  _rE(pthread_attr_setscope(&threadAttr, PTHREAD_SCOPE_SYSTEM),          "sweatShop::run()--  Failed to configure pthreads (set scope): %s.\n");
  _rE(pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_JOINABLE), "sweatShop::run()--  Failed to configure pthreads (joinable): %s.\n");
  _rE(pthread_create(&threadIDloader, &threadAttr, _sload, this),        "sweatShop::run()--  Failed to launch loader thread: %s.\n");

  while (!_writerP && !_workerP && !_loaderP)   //  Wait for the loader to actually load data,
    nanosleep(&naptime, nullptr);               //  otherwise all the workers immediately quit.

  //  Start the statistics and writer

  _rE(pthread_create(&threadIDstats,  &threadAttr, _sstat, this), "sweatShop::run()--  Failed to launch status thread: %s.\n");
  _rE(pthread_create(&threadIDwriter, &threadAttr, _swrit, this), "sweatShop::run()--  Failed to launch writer thread: %s.\n");

  //  And some labor

  for (uint32 i=0; i<_numberOfWorkers; i++)
    _rE(pthread_create(&_workerData[i].threadID, &threadAttr, _swork, _workerData + i), "sweatShop::run()--  Failed to launch worker thread " F_U32 ": %s.\n", i);

  //  Now sit back and relax.

  _rE(pthread_join(threadIDloader, nullptr), "sweatShop::run()--  Failed to join loader thread: %s.\n");
  _rE(pthread_join(threadIDwriter, nullptr), "sweatShop::run()--  Failed to join writer thread: %s.\n");
  _rE(pthread_join(threadIDstats,  nullptr), "sweatShop::run()--  Failed to join status thread: %s.\n");

  for (uint32 i=0; i<_numberOfWorkers; i++)
    _rE(pthread_join(_workerData[i].threadID, nullptr), "sweatShop::run()--  Failed to join worker thread " F_U32 ": %s.\n", i);

  //  Cleanup.

  pthread_attr_destroy(&threadAttr);

  for (uint32 i=0; i<_numberOfWorkers; i++)
    delete [] _workerData[i].workerQueue;

  delete _loaderP;

  _loaderP = _workerP = _writerP = nullptr;
}
