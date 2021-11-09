
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

#include "meryl.H"
#include "kmers.H"
#include "sequence.H"

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
  ~merylInput();

  char const *inputType(void) {
    if (_operation)   return("meryl-operation");
    if (_stream)      return("meryl-database");
    if (_sequence)    return("sequence-file");
    if (_store)       return("canu-seqStore");

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

  bool   isMultiSet(void)         { return(_isMultiSet);        };

  bool   isCompressedFile(void)   { return((_sequence != NULL) &&
                                           (_sequence->isCompressed())); };

  merylOperation           *_operation;
  merylFileReader          *_stream;
  dnaSeqFile               *_sequence;
  sqStore                  *_store;

  char                      _name[FILENAME_MAX+1];
  bool                      _isMultiSet;

  //  For _operation and _stream, a copy of the 'active' kmer

  kmer                      _kmer;
  kmvalu                    _value;
  bool                      _valid;

  //  For _store, the position we're at in the store.

  uint32                    _sqBgn;
  uint32                    _sqEnd;

  sqRead                   *_read;
  uint32                    _readID;
  uint32                    _readPos;

  //  For homopoly compression, a flag, and the last base read.

  bool                      _homopolyCompress;
  char                      _lastByte;
};

#endif  //  MERYLINPUT_H
