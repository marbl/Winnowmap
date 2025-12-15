
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

#ifndef MERYLUTIL_SEQUENCE_V1_DNASEQFILE_H
#define MERYLUTIL_SEQUENCE_V1_DNASEQFILE_H

#include "types.H"
#include "arrays.H"
#include "files.H"

#include "htslib/hts/sam.h"

#include "dnaSeq-v1.H"

//
//  An interface to FASTA, FASTQ, SAM, BAM and CRAM files.
//
//  Open a sequence file with:
//    dnaSeqFile *f = openSequenceFile(name, wantIndex);
//
//  If the file cannot be loaded, nullptr is returned.
//
//  Without an index, numberOfSequences(), findSequence() and
//  sequenceLength() do not work well or at all.  SAM, BAM and CRAM files do
//  NOT support indexes.
//
//  You can manually load, create or remove (from memory) an index with
//    loadIndex()    - load index into memory, option to create if needed
//    unloadIndex()  - remove index from memory
//    createIndex()  - build on-disk index
//    destroyIndex() - remove on-disk index
//
//  close() will close the file; useful only for internal operations.
//  Users should just delete the object.
//
//  reopen() behaves as expected; it is effectively a rewind() operation, but
//  explicitly closed and opens the file again.  If an index was loaded
//  before, it is loaded again.
//
//  Operations to get sequence from a dnaSeqFile:
//
//  numberOfSequences()
//   - if an index exists, will return the number of sequences in the file,
//     otherwise, returns 0.
//
//  findSequence()
//   - returns true if the specified sequence is found in the file and will
//     leave the file positioned such that the next loadSequence() will load
//     that sequence.
//   - If an index exists, the index will be searched and the sequence will
//     be returned regardless of where it is in the file.
//   - If no index exists, the file will be searched forward until the
//     sequence is found or the file ends.  It is not possible to move
//     'backward' in the file in this case.
//
//  sequenceLength()
//   - if an index exists, will return the length of sequence index i,
//     otherwise, returns uint64max.
//
//  isFASTA() and isFASTQ()
//   - return true if the last sequence loaded came from a FASTA or FASTQ
//     source, respectively.  If no sequence has been loaded yet, both
//     functions will return false.
//   - only generally useful internally, for outputting a sequence in the
//     same format it was loaded from disk.
//
//  loadSequence()
//   - will read the next sequence from the file.  Returns false if the end
//     of file is encountered, true otherwise.  In particular, a sequence of
//     length zero will return true.
//
//  loadBases()
//   - will return a chunk of sequence from the file, up to 'maxLength' bases
//     or the end of the current sequence.
//   - Returns false only if EOF is encountered.
//   - seqLength will have the length of the sequence returned.  This can be zero.
//   - endOfSequence will be true if the end of the sequence was encountered.
//   - The returned sequence is NOT NUL terminated.
//

namespace merylutil::inline sequence::inline v1 {

class dnaSeqFile {
public:
  dnaSeqFile()  {                }
  virtual
  ~dnaSeqFile() { delete [] _fn; }

  void            filename(char const *f) { _fn = duplicateString(f); }
  char const     *filename(void)          { return _fn;   }

  //  Compressed and streaming inputs are not indexable.
  //  Streaming inputs are not reopenable.
  bool            isIndexable(void)     { return _indexable;  }
  bool            isReopenable(void)    { return _reopenable; }
  bool            isCompressed(void)    { return _compressed; }

  //rtual bool    open(char const *filename, bool indexed=false) = 0;
  virtual bool    close(void) = 0;
  virtual bool    reopen(void) = 0;

  virtual bool    loadIndex(bool create=false) = 0;
  virtual void    unloadIndex(void) = 0;
  virtual void    createIndex(void) = 0;
  virtual void    destroyIndex(void) = 0;

  //  Operations requiring an index.

  uint64          numberOfSequences(void)  { return _numSeqs; }

  virtual bool    findSequence(uint64 i) = 0;
  virtual uint64  sequenceLength(uint64 i) = 0;

  //  Operations returning sequence

  virtual bool    loadSequence(char   *&name, uint32 &nameMax,
                               char   *&seq,
                               uint8  *&qlt,  uint64 &seqMax, uint64 &seqLen, uint32 &errorCode) = 0;

  bool            loadSequence(dnaSeq &seq) {
    if (loadSequence(seq._headr, seq._headrMax,
                     seq._seq,
                     seq._qlt, seq._seqMax, seq._seqLen, seq._error)) {
      seq.findNameAndFlags();
      return true;
    }
    else
      return false;
  };

  bool            isFASTA(void)      { return(_isFASTA);  }   //  The last sequence returned
  bool            isFASTQ(void)      { return(_isFASTQ);  }   //  is FASTA/Q, and has index
  uint64          seqIdx(void)       { return(_seqIdx-1); }   //  seqIdx().

  virtual bool    loadBases(char *seq, uint64 maxLength, uint64 &seqLength, bool &endOfSequence) = 0;

private:
  char      *_fn           = nullptr;
protected:
  bool       _indexable    = false;   //  True if we can randomly access the file.
  bool       _reopenable   = false;   //  True if we can reopen the file (e.g., not stdin).
  bool       _compressed   = false;   //  True if a second process is running to decompress the file.

  uint64     _numSeqs      = 0;

  bool       _isFASTA      = false;
  bool       _isFASTQ      = false;
  uint64     _seqIdx       = 0;
};



class bufSeqFile : public dnaSeqFile {
public:
  bufSeqFile()  {          }
  ~bufSeqFile() { close(); }

  bool        open(char const *filename, bool indexed=false);
  bool        close(void);
  bool        reopen(void);

  bool        loadIndex(bool create=true);
  void        unloadIndex(void);
  void        createIndex(void);
  void        destroyIndex(void);

  bool        findSequence(uint64 i);
  uint64      sequenceLength(uint64 i);

private:
  bool
  loadFASTA(char  *&name, uint32 &nameMax,
            char  *&seq,
            uint8 *&qlt,  uint64 &seqMax, uint64 &seqLen, uint64 &qltLen);

  bool
  loadFASTQ(char  *&name, uint32 &nameMax,
            char  *&seq,
            uint8 *&qlt,  uint64 &seqMax, uint64 &seqLen, uint64 &qltLen);

public:
  bool        loadSequence(char   *&name, uint32 &nameMax,
                           char   *&seq,
                           uint8  *&qlt,  uint64 &seqMax, uint64 &seqLen, uint32 &errorCode);
  //bool        loadSequence(dnaSeq &seq)

  bool        loadBases(char *seq, uint64 maxLength, uint64 &seqLength, bool &endOfSequence);



private:
  compressedFileReader  *_file         = nullptr;
  readBuffer            *_buffer       = nullptr;

  struct dsfIdxEnt {            //  Offset of the first byte in the record:
    uint64   _fileOffset;       //  '>' for FASTA, '@' for fastq.
    uint64   _sequenceLength;   //
  };

  dsfIdxEnt *_index        = nullptr;
  uint64     _indexLen     = 0;
  uint64     _indexMax     = 0;
};



class htsSeqFile : public dnaSeqFile {
public:
  htsSeqFile()  {          }
  ~htsSeqFile() { close(); }

  bool        open(char const *filename, bool indexed=false);
  bool        close(void);
  bool        reopen(void);

  bool        loadIndex(bool create=true);
  void        unloadIndex(void);
  void        createIndex(void);
  void        destroyIndex(void);

  bool        findSequence(uint64 i);
  uint64      sequenceLength(uint64 i);

  bool        loadSequence(char   *&name, uint32 &nameMax,
                           char   *&seq,
                           uint8  *&qlt,  uint64 &seqMax, uint64 &seqLen, uint32 &errorCode);
  //bool        loadSequence(dnaSeq &seq);

  bool        loadBases(char *seq, uint64 maxLength, uint64 &seqLength, bool &endOfSequence);

private:
  htsFile    *_hts    = nullptr;
  sam_hdr_t  *_htshdr = nullptr;
  bam1_t     *_htsbam = nullptr;
  uint64      _htspos = uint64max;
};


inline
dnaSeqFile *
openSequenceFile(char const *filename, bool index=false) {

  bufSeqFile  *bs = new bufSeqFile();
  if (bs->open(filename, index))
    return bs;
  else
    delete bs;

  htsSeqFile  *hs = new htsSeqFile();
  if (hs->open(filename, index))
    return hs;
  else
    delete hs;

  fprintf(stderr, "Don't know how to open sequence file '%s'.\n", filename);

  return nullptr;
}


}  //  namespace merylutil::sequence::v1

#endif  //  MERYLUTIL_SEQUENCE_V1_DNASEQFILE_H
