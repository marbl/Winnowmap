
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

#include "meryl.H"
#include "strings.H"
#include "system.H"
#include "sweatShop.H"

#include <atomic>


class mcGlobalData {
public:
  mcGlobalData(std::vector<merylInput *> &inputs,
               merylOp                    op,
               uint64                     nPrefix,
               uint32                     wData,
               kmdata                     wDataMask,
               uint64                     maxMemory,
               uint32                     maxThreads,
               uint64                     bufferSize,
               merylFileWriter           *output) : _inputs(inputs) {
    _operation      = op;
    _nPrefix        = nPrefix;
    _wData          = wData;
    _wDataMask      = wDataMask;

    _dumping        = false;

    _lock           = new std::atomic_flag [_nPrefix];
    _data           = new merylCountArray  [_nPrefix];
    _output         = output;
    _writer         = output->getBlockWriter();

    _maxMemory      = maxMemory;
    _memBase        = getProcessSize();
    _memUsed        = _memBase;
    _memReported    = 0;

    _maxThreads     = maxThreads;
    _loadThreads    = 1;

    _bufferSize     = bufferSize;

    _kmersAdded     = 0;
    _kmersAddedMax  = 0;

    _inputPos       = 0;

    for (uint32 ii=0; ii<65; ii++)
      _lastBuffer[ii] = 0;

    for (uint32 pp=0; pp<_nPrefix; pp++) {      //  Initialize each bucket.
      _lock[pp].clear();
      _memUsed += _data[pp].initialize(pp, wData);
    }
  };

  ~mcGlobalData() {
    delete [] _lock;
    delete [] _data;
    delete [] _writer;
  };

  merylOp                     _operation;        //  Parameters.
  uint64                      _nPrefix;
  uint32                      _wData;
  kmdata                      _wDataMask;

  bool                        _dumping;

  std::atomic_flag           *_lock;
  merylCountArray            *_data;             //  Data for counting.
  merylFileWriter            *_output;
  merylBlockWriter           *_writer;           //  Data for writing.

  uint64                      _maxMemory;        //  Maximum memory we can use.
  uint64                      _memBase;          //  Overhead memory.
  uint64                      _memUsed;          //  Sum of actual memory used.
  uint64                      _memReported;      //  Memory usage at last report.

  uint32                      _maxThreads;       //  The max number of CPUs we can use.
  uint32                      _loadThreads;      //  The number of CPUs used for reading input.

  uint64                      _bufferSize;       //  Maximum size of a computation input buffer.

  uint64                      _kmersAdded;       //  Number of kmers added; boring statistics for the user.
  uint64                      _kmersAddedMax;    //  Max kmers in any single merylCountArray; not boring.
  
  uint32                      _inputPos;         //  Input files.
  std::vector<merylInput *>  &_inputs;

  char                        _lastBuffer[65];   //  Wrap-around from the last buffer.
};



class mcComputation {
public:
  mcComputation(uint64 bufmax) {
    _bufferMax     = bufmax;
    _buffer        = new char [_bufferMax];
  };

  ~mcComputation() {
    delete [] _buffer;
  };

  uint64        _bufferMax  = 0;       //  Input data
  uint64        _bufferLen  = 0;
  char         *_buffer     = nullptr;

  kmerIterator  _kiter;                //  Sequence to kmer conversion

  uint64        _memUsed       = 0;    //  Output statistics on kmers added to
  uint64        _kmersAdded    = 0;    //  the merylCountArray but this block.
  uint64        _kmersAddedMax = 0;
};




void *
loadBases(void *G) {
  mcGlobalData     *g  = (mcGlobalData  *)G;
  mcComputation    *s  = new mcComputation(g->_bufferSize);
  uint32            kl = kmerTiny::merSize() - 1;

  //  Copy the end of the last block into our buffer.

  assert(s->_bufferLen == 0);
  assert(s->_bufferMax > kl);

  if (g->_lastBuffer[0] != 0) {
    memcpy(s->_buffer, g->_lastBuffer, sizeof(char) * kl);

    s->_bufferLen += kl;

    g->_lastBuffer[0] = 0;
  }

  //  If no more inputs, we're done.

  if (g->_inputPos >= g->_inputs.size())
    return(NULL);

  //  Update the number of threads used for loading.  If the input is
  //  compressed, reserve 2 threads, otherwise reserve 1.

  if (g->_inputs[g->_inputPos]->isCompressedFile())
    g->_loadThreads = 2;
  else
    g->_loadThreads = 1;

  //  Try to load bases.  Keep loading until the buffer is filled
  //  or we exhaust the file.

  while (1) {
    uint64  bMax     = s->_bufferMax - s->_bufferLen;
    uint64  bLen     = 0;
    bool    endOfSeq = false;

    if (bMax < 512)   //  If the buffer is full enough,
      break;          //  stop loading.

    //  Load bases, but reserve 2 characters in the buffer for a
    //  sequence terminating ','.

    bool success = g->_inputs[g->_inputPos]->loadBases(s->_buffer + s->_bufferLen,
                                                       bMax - 2,
                                                       bLen, endOfSeq);

    //  If no bases loaded, we've exhausted the file.  Close it,
    //  and move to the next one.  Bail on loading more bases;
    //  let the test above figure out if we're all done on
    //  the next call to loadBases().

    if (success == false) {
      assert(bLen == 0);

      s->_buffer[s->_bufferLen++] = '.';   //  Insert a mer-breaker, just to be safe.

      delete g->_inputs[g->_inputPos]->_sequence;
      g->_inputs[g->_inputPos]->_sequence = NULL;

      g->_inputPos++;

      break;
    }

    //  Account for whatever we just loaded.

    s->_bufferLen += bLen;

    assert(s->_bufferLen+1 <= s->_bufferMax);

    //  If we're at the end of a sequence, append a mer-breaker.

    if (endOfSeq == true)
      s->_buffer[s->_bufferLen++] = '.';
  }

  //  With bases loaded, we need to save the last few bases for the next buffer,
  //  and tell the kmerIterator about the bases we loaded.

  if (s->_buffer[s->_bufferLen-1] != '.')
    memcpy(g->_lastBuffer, s->_buffer + s->_bufferLen - kl, sizeof(char) * kl);

  //  Now just tell the iterator about the buffer.

  s->_kiter.addSequence(s->_buffer, s->_bufferLen);

  //  That's it.

  return(s);
}



void
insertKmers(void *G, void *T, void *S) {
  mcGlobalData     *g = (mcGlobalData  *)G;
  mcComputation    *s = (mcComputation *)S;

  while (s->_kiter.nextMer()) {
    bool    useF = (g->_operation == opCountForward);
    kmdata  pp   = 0;
    kmdata  mm   = 0;

    if (g->_operation == opCount)
      useF = (s->_kiter.fmer() < s->_kiter.rmer());

    if (useF == true) {
      pp = (kmdata)s->_kiter.fmer() >> g->_wData;
      mm = (kmdata)s->_kiter.fmer()  & g->_wDataMask;
      //fprintf(stderr, "useF F=%s R=%s ms=%u pp %llu mm %llu\n", s->_kiter.fmer().toString(fstr), s->_kiter.rmer().toString(rstr), s->_kiter.fmer().merSize(), pp, mm);
    }

    else {
      pp = (kmdata)s->_kiter.rmer() >> g->_wData;
      mm = (kmdata)s->_kiter.rmer()  & g->_wDataMask;
      //fprintf(stderr, "useR F=%s R=%s ms=%u pp %llu mm %llu\n", s->_kiter.fmer().toString(fstr), s->_kiter.rmer().toString(rstr), s->_kiter.rmer().merSize(), pp, mm);
    }

    assert(pp < g->_nPrefix);

    //  If we're dumping data, stop immediately and sleep until dumping is
    //  finished.

    while (g->_dumping == true)
      usleep(1000);

    //  We need exclusive access to this specific merylCountArray, so busy
    //  wait on a lock until we get it.

    while (g->_lock[pp].test_and_set(std::memory_order_relaxed) == true)
      ;

    s->_memUsed        += g->_data[pp].add(mm);
    s->_kmersAdded     += 1;
    s->_kmersAddedMax   = std::max(s->_kmersAddedMax, g->_data[pp].numKmers());

    g->_lock[pp].clear(std::memory_order_relaxed);
  }
}




void
writeBatch(void *G, void *S) {
  mcGlobalData     *g = (mcGlobalData  *)G;
  mcComputation    *s = (mcComputation *)S;

  //  Udpate memory used and kmers added.  There's only one writer thread,
  //  so this is thread safe!

  g->_memUsed       += s->_memUsed;
  g->_kmersAdded    += s->_kmersAdded;
  g->_kmersAddedMax  = std::max(s->_kmersAddedMax, g->_kmersAddedMax);

  //  Free the input buffer.  All the data is loaded into merylCountArrays,
  //  and all we needed to get from this is the stats above.

  delete s;

  //  Estimate, poorly, how much memory we'll need to sort the arrays.  It's
  //  a poor estimate because we'll never have all threads sorting the
  //  maximum number of kmers at the same time, but it's a safe poor
  //  estimate.

  uint64  sortMem = g->_maxThreads * g->_kmersAddedMax * sizeof(kmdata);

  //  Write a log every 128 MB of memory growth.

  if (g->_memUsed + sortMem - g->_memReported > (uint64)128 * 1024 * 1024) {
    g->_memReported = g->_memUsed + sortMem;

    fprintf(stderr, "Used %3.3f GB / %3.3f GB to store %12lu kmers; need %3.3f GB to sort %12lu kmers\n",
            g->_memUsed   / 1024.0 / 1024.0 / 1024.0,
            g->_maxMemory / 1024.0 / 1024.0 / 1024.0,
            g->_kmersAdded,
            sortMem / 1024.0 / 1024.0 / 1024.0, g->_kmersAddedMax);
  }

  //  If we haven't hit the memory limit yet, just return.

  if (g->_memUsed + sortMem < g->_maxMemory)
    return;

  //  Tell all the threads to pause, then grab all the locks to ensure nobody
  //  is still adding kmers to a merylCountArray.

  g->_dumping = true;

  for (uint32 pp=0; pp<g->_nPrefix; pp++)
    while (g->_lock[pp].test_and_set(std::memory_order_relaxed) == true)
      ;

  //  Write data!  For reasons I don't understand, we need to reset the max
  //  number of threads to use.  Something is resetting it to the number of
  //  CPUs on the machine.
  //
  //  Since we still have a sequence loader around, we need to leave threads
  //  for it.

  uint32  wThreads = (g->_maxThreads > g->_loadThreads) ? (g->_maxThreads - g->_loadThreads) : 1;
  uint32  lThreads =                   g->_loadThreads;

  fprintf(stderr, "Memory full.  Writing results to '%s', using %u thread%s (%u thread%s still doing input).\n",
          g->_output->filename(),
          wThreads, (wThreads == 1) ? "" : "s",
          lThreads, (lThreads == 1) ? "" : "s");
  fprintf(stderr, "\n");

  setNumThreads(wThreads);

#pragma omp parallel for schedule(dynamic, 1)
  for (uint32 ff=0; ff<g->_output->numberOfFiles(); ff++) {
    for (uint64 pp=g->_output->firstPrefixInFile(ff); pp <= g->_output->lastPrefixInFile(ff); pp++) {
      g->_data[pp].countKmers();                   //  Convert the list of kmers into a list of (kmer, count).
      g->_data[pp].dumpCountedKmers(g->_writer);   //  Write that list to disk.
      g->_data[pp].removeCountedKmers();           //  And remove the in-core data.
    }
  }

  g->_writer->finishBatch();

  //  Reset accounting.

  g->_memUsed    = g->_memBase;

  for (uint32 pp=0; pp<g->_nPrefix; pp++)
    g->_memUsed += g->_data[pp].usedSize();

  g->_kmersAdded    = 0;
  g->_kmersAddedMax = 0;

  //  Signal that threads can proceeed.

  for (uint32 pp=0; pp<g->_nPrefix; pp++)
    g->_lock[pp].clear(std::memory_order_relaxed);

  g->_dumping = false;
}



void
merylOperation::countThreads(uint32  wPrefix,
                             uint64  nPrefix,
                             uint32  wData,
                             kmdata  wDataMask) {

  //  If we're only configuring, stop now.

  if (_onlyConfig)
    return;

  //  Configure the writer for the prefix bits we're counting with.
  //
  //  We split the kmer into wPrefix and wData (bits) pieces.
  //  The prefix is used by the filewriter to decide which file to write to, by simply
  //  shifting it to the right to keep the correct number of bits in the file.
  //
  //           kmer -- [ wPrefix (18) = prefixSize               | wData (36) ]
  //           file -- [ numFileBits  | prefixSize - numFileBits ]

  _outputO->initialize(wPrefix);

  //  Initialize the counter.
  //
  //  Tell it to use _maxMemory, but carve out space for the input buffers.
  //  At 2MB each, 16 per thread, and 16 threads, that's 512 MB.  Not huge,
  //  but a big chunk of our (expected 16 or so GB total).  The extra buffers
  //  are generally filled when a batch is dumped to disk.

  uint64  inputBufferSize = 2 * 1024 * 1024;

  mcGlobalData  *g = new mcGlobalData(_inputs,
                                      _operation,
                                      nPrefix,
                                      wData,
                                      wDataMask,
                                      _maxMemory - inputBufferSize * 4 * _maxThreads,
                                      _maxThreads,
                                      inputBufferSize,
                                      _outputO);

  //  Set up a sweatShop and run it.  We'll reserve one thread for input, one
  //  for gzip and use the remaining for counting -- unless there are no
  //  remaining, then we'll just use one.

  sweatShop    *ss = new sweatShop(loadBases, insertKmers, writeBatch);

  uint32 nw = (_maxThreads > 2) ? (_maxThreads - 2) : 1;

  ss->setLoaderBatchSize(1);            //  Load this many things before appending to input list
  ss->setLoaderQueueSize(nw * 16);      //  Allow this many things on the input list before stalling the input
  ss->setWriterQueueSize(nw);           //  Allow this many things on the output list before stalling the compute
  ss->setNumberOfWorkers(nw);           //  Use this many worker CPUs; leave one for input and one for gzip.

  ss->run(g, false);

  delete ss;

  //  All data loaded.  Write the output.  Reset threads before starting (see
  //  above) to the maximum possible since there is no loader threads around
  //  anymore.

  fprintf(stderr, "\n");
  fprintf(stderr, "Input complete.  Writing results to '%s', using %u thread%s.\n",
          _outputO->filename(), _maxThreads, (_maxThreads == 1) ? "" : "s");

  setNumThreads(_maxThreads);

#pragma omp parallel for schedule(dynamic, 1)
  for (uint32 ff=0; ff<_outputO->numberOfFiles(); ff++) {
    for (uint64 pp=_outputO->firstPrefixInFile(ff); pp <= _outputO->lastPrefixInFile(ff); pp++) {
      g->_data[pp].countKmers();                   //  Convert the list of kmers into a list of (kmer, count).
      g->_data[pp].dumpCountedKmers(g->_writer);   //  Write that list to disk.
      g->_data[pp].removeCountedKmers();           //  And remove the in-core data.
    }
  }

  //  Merge any iterations into a single file, or just rename
  //  the single file to the final name.

  g->_writer->finish();

  delete g->_writer;   //  Explicitly delete the writer.
  g->_writer = NULL;

  //  Cleanup.
  delete g;

  fprintf(stderr, "\n");
  fprintf(stderr, "Finished counting.\n");
}
