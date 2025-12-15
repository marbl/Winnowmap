
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

#ifndef MERYLINPUT_H
#define MERYLINPUT_H

#include "sequence.H"
#include "kmers.H"

using namespace merylutil;
using namespace merylutil::kmers::v1;

#ifdef CANU                   //  If we're building under
#include "sqStore.H"          //  Canu, include support
#else                         //  to read the seqStore
typedef  uint32  sqStore;     //  directly, otherwise,
typedef  uint32  sqRead;      //  add fake types to reduce
#endif                        //  the number of ifdef blocks.

class merylOperation;


class merylInput {
public:
  merylInput(merylOperation *o);
  merylInput(const char *n, merylFileReader *s, uint32 fileNum=UINT32_MAX);
  merylInput(const char *n, dnaSeqFile *f, bool doCompression);
  merylInput(const char *n, sqStore *s, uint32 segment, uint32 segmentMax);
  merylInput(const char *n);
  ~merylInput();

  char const *inputType(void) {
    if (_operation)   return("meryl-operation");
    if (_stream)      return("meryl-database");
    if (_sequence)    return("sequence-file");
    if (_store)       return("canu-seqStore");
    if (_filename)    return("text-file");

    return("invalid-input");
  };

  void   initialize(void);
  void   nextMer(void);

private:
  bool   loadBasesFromCanu(char    *seq,
                           uint64   maxLength,
                           uint64  &seqLength,
                           bool    &endOfSequence);

public:
  bool   loadBases(char    *seq,
                   uint64   maxLength,
                   uint64  &seqLength,
                   bool    &endOfSequence);

  bool   isFromOperation(void)    { return(_operation != NULL); };
  bool   isFromDatabase(void)     { return(_stream    != NULL); };
  bool   isFromSequence(void)     { return(_sequence  != NULL); };
  bool   isFromStore(void)        { return(_store     != NULL); };
  bool   isFromFile(void)         { return(_filename  != NULL); };

  bool   isMultiSet(void)         { return(_isMultiSet);        };

  bool   isCompressed(void)       { return((_sequence != NULL) &&
                                           (_sequence->isCompressed())); };

public:
  char                      _name[FILENAME_MAX+1] = {0};

  merylOperation           *_operation = nullptr;
  merylFileReader          *_stream    = nullptr;
  dnaSeqFile               *_sequence  = nullptr;
  sqStore                  *_store     = nullptr;
  const char               *_filename  = nullptr;

  bool                      _isMultiSet = false;

  //  For _operation and _stream, a copy of the 'active' kmer

  kmer                      _kmer;
  kmvalu                    _value  = 0;
  bool                      _valid  = false;

  //  For _store, the position we're at in the store.

  uint32                    _sqBgn = 0;
  uint32                    _sqEnd = 0;

  sqRead                   *_read    = nullptr;
  uint32                    _readID  = 0;
  uint32                    _readPos = uint32max;

  //  For homopoly compression, a flag, and the last base read.

  bool                      _homopolyCompress = false;
  char                      _lastByte = 0;
};

#endif  //  MERYLINPUT_H
