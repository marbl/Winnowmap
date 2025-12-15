

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



void
merylOpCounting::countSequential(std::vector<merylInput *> &inputs,
                                 uint64                     allowedMemory,
                                 uint32                     allowedThreads,
                                 merylFileWriter           *output) {

  //  If we're only configuring, stop now.

  if (_onlyConfig)
    return;

  //  Configure the writer for the prefix bits we're counting with.
  //
  //  We split the kmer into _wPrefix and _wSuffix (bits) pieces.
  //  The prefix is used by the filewriter to decide which file to write to, by simply
  //  shifting it to the right to keep the correct number of bits in the file.
  //
  //           kmer -- [ _wPrefix (18) = prefixSize               | _wSuffix (36) ]
  //           file -- [ numFileBits  | prefixSize - numFileBits ]

  output->initialize(_wPrefix);

  merylBlockWriter  *_writer = output->getBlockWriter();

  //  Allocate buckets.  The buckets don't allocate space for mers until they're added,
  //  and allocate space for these mers in blocks of 8192 * 64 bits.
  //
  //  Need someway of balancing the number of prefixes we have and the size of each
  //  initial allocation.

  merylCountArray  *data = new merylCountArray [_nPrefix];

  //  Load bases, count!

  uint64          bufferMax  = 1024 * 1024;
  uint64          bufferLen  = 0;
  char           *buffer     = new char [bufferMax];
  bool            endOfSeq   = false;

  memset(buffer, 0, sizeof(char) * bufferMax);

  //char            fstr[65];
  //char            rstr[65];

  uint64          memBase     = getProcessSize();   //  Overhead memory.
  uint64          memUsed     = 0;                  //  Sum of actual memory used.
  uint64          memReported = 0;                  //  Memory usage at last report.

  memUsed = memBase;

  for (uint32 pp=0; pp<_nPrefix; pp++)
    memUsed += data[pp].initialize(pp, _wSuffix);

  uint64          kmersAdded  = 0;

  kmerIterator    kiter;

  for (uint32 ii=0; ii<inputs.size(); ii++) {
    fprintf(stderr, "Loading kmers from '%s' into buckets.\n", inputs[ii]->inputName());

    while (inputs[ii]->loadBases(buffer, bufferMax, bufferLen, endOfSeq)) {
      if (bufferLen == 0)
        continue;

      //fprintf(stderr, "read " F_U64 " bases from '%s'\n", bufferLen, inputs[ii]->_name);

      kiter.addSequence(buffer, bufferLen);

      while (kiter.nextMer()) {
        bool    useF = _countForward;
        kmdata  pp   = 0;
        kmdata  mm   = 0;

        if (_countCanonical == true)
          useF = (kiter.fmer() < kiter.rmer());

        if (useF == true) {
          pp = (kmdata)kiter.fmer() >> _wSuffix;
          mm = (kmdata)kiter.fmer()  & _wSuffixMask;
          //fprintf(stderr, "useF F=%s R=%s ms=%u pp %lu mm %lu\n", kiter.fmer().toString(fstr), kiter.rmer().toString(rstr), kiter.fmer().merSize(), pp, mm);
        }

        else {
          pp = (kmdata)kiter.rmer() >> _wSuffix;
          mm = (kmdata)kiter.rmer()  & _wSuffixMask;
          //fprintf(stderr, "useR F=%s R=%s ms=%u pp %lu mm %lu\n", kiter.fmer().toString(fstr), kiter.rmer().toString(rstr), kiter.rmer().merSize(), pp, mm);
        }

        assert(pp < _nPrefix);

        memUsed += data[pp].add(mm);

        kmersAdded++;
      }

      if (endOfSeq)      //  If the end of the sequence, clear
        kiter.reset();   //  the running kmer.

      //  Report that we're actually doing something.

      if (memUsed - memReported > (uint64)128 * 1024 * 1024) {
        memReported = memUsed;

        fprintf(stderr, "Used %3.3f GB out of %3.3f GB to store %12lu kmers.\n",
                memUsed       / 1024.0 / 1024.0 / 1024.0,
                allowedMemory / 1024.0 / 1024.0 / 1024.0,
                kmersAdded);
      }

      //  If we're out of space, process the data and dump.

      if (memUsed > allowedMemory) {
        fprintf(stderr, "Memory full.  Writing results to '%s', using " F_S32 " threads.\n",
                output->filename(), getMaxThreadsAllowed());
        fprintf(stderr, "\n");

#pragma omp parallel for schedule(dynamic, 1)
        for (uint32 ff=0; ff<output->numberOfFiles(); ff++) {
          //fprintf(stderr, "thread %2u writes file %2u with prefixes 0x%016lx to 0x%016lx\n",
          //        omp_get_thread_num(), ff, output->firstPrefixInFile(ff), output->lastPrefixInFile(ff));

          for (uint64 pp=output->firstPrefixInFile(ff); pp <= output->lastPrefixInFile(ff); pp++) {
            data[pp].countKmers();                            //  Convert the list of kmers into a list of (kmer, count).
            data[pp].dumpCountedKmers(_writer, _lConstant);   //  Write that list to disk.
            data[pp].removeCountedKmers();                    //  And remove the in-core data.
          }
        }

        _writer->finishBatch();

        kmersAdded = 0;

        memUsed = memBase;                        //  Reinitialize or memory used.
        for (uint32 pp=0; pp<_nPrefix; pp++)
          memUsed += data[pp].usedSize();
      }

    }

    //  Would like some kind of report here on the kmers loaded from this file.

    delete inputs[ii]->_sequence;
    inputs[ii]->_sequence = NULL;
  }

  //  Finished loading kmers.  Free up some space.

  //delete [] kmers;
  delete [] buffer;

  //  Sort, dump and erase each block.
  //
  //  A minor complication is that within each output file, the blocks must be in order.

  fprintf(stderr, "\n");
  fprintf(stderr, "Writing results to '%s', using " F_S32 " threads.\n",
          output->filename(), getMaxThreadsAllowed());

  //for (uint64 pp=0; pp<_nPrefix; pp++)
  //  fprintf(stderr, "Prefix 0x%016lx writes to file %u\n", pp, output->fileNumber(pp));

#pragma omp parallel for schedule(dynamic, 1)
  for (uint32 ff=0; ff<output->numberOfFiles(); ff++) {
    //fprintf(stderr, "thread %2u writes file %2u with prefixes 0x%016lx to 0x%016lx\n",
    //        omp_get_thread_num(), ff, output->firstPrefixInFile(ff), output->lastPrefixInFile(ff));

    for (uint64 pp=output->firstPrefixInFile(ff); pp <= output->lastPrefixInFile(ff); pp++) {
      data[pp].countKmers();                            //  Convert the list of kmers into a list of (kmer, count).
      data[pp].dumpCountedKmers(_writer, _lConstant);   //  Write that list to disk.
      data[pp].removeCountedKmers();                    //  And remove the in-core data.
    }
  }

  //  Merge any iterations into a single file, or just rename
  //  the single file to the final name.

  _writer->finish();

  delete _writer;
  _writer = NULL;

  //  Cleanup.

  delete [] data;

  fprintf(stderr, "\n");
  fprintf(stderr, "Finished counting.\n");
}
