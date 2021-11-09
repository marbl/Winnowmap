
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

#ifndef MERYL_UTIL_KMER_FILES_H
#define MERYL_UTIL_KMER_FILES_H

#ifndef MERYL_UTIL_KMER_H
#error "include kmers.H, not this."
#endif



void
dumpMerylDataFile(char *name);



//  Functions to constrct data file names and open them for reading or
//  writing.

char *
constructBlockName(char   *nameprefix,
                   uint64  outIndex,
                   uint32  numFiles,
                   uint32  iteration,
                   bool    isIndex);

FILE *
openOutputBlock(char   *nameprefix,
                uint64  fileIndex,
                uint32  numFiles,
                uint32  iteration=0);

FILE *
openInputBlock(char   *nameprefix,
               uint64  fileIndex,
               uint32  numFiles,
               uint32  iteration=0);


//  Read a block of kmer data from disk, and decode it into a list of kmers,
//  counts and (eventually) colors.

class merylFileBlockReader {
public:
  merylFileBlockReader() {
  }
  ~merylFileBlockReader() {
    delete [] _suffixes;
    delete [] _values;
  }

  bool      loadBlock(FILE *inFile, uint32 activeFile, uint32 activeIteration=0);

  void      decodeBlock(void);                               //  to our own storage
  void      decodeBlock(kmdata *suffixes, kmvalu *values);   //  to external storage

  kmpref    prefix(void)   { return(_blockPrefix); };        //  kmer prefix of this block
  uint64    nKmers(void)   { return(_nKmers);      };        //  number of kmers in this block

  kmdata   *suffixes(void) { return(_suffixes); };           //  direct access to decoded data
  kmvalu   *values(void)   { return(_values);   };

private:
  stuffedBits   _data;

  kmpref        _blockPrefix = 0;   //  The prefix of all kmers in this block
  uint64        _nKmers      = 0;   //  The number of kmers in this block
  uint64        _nKmersMax   = 0;   //  The number of kmers we've allocated space for in _suffixes and _values

  uint32        _kCode       = 0;   //  Encoding type of kmer, then 128 bits of parameters
  uint32        _unaryBits   = 0;   //    bits in the unary prefix  (of the kmer suffix)
  uint32        _binaryBits  = 0;   //    bits in the binary suffix (of the kmer suffix)
  uint64        _k1          = 0;   //    unused

  uint32        _cCode       = 0;   //  Encoding type of the values, then 128 bits of parameters
  uint64        _c1          = 0;   //    unused
  uint64        _c2          = 0;   //    unused

  kmdata       *_suffixes    = nullptr;   //  Decoded suffixes and values.
  kmvalu       *_values      = nullptr;   //
};


//  An index to the binary encoded kmer data.  Provides:
//    kmer prefix for each block
//    starting position of the block in the file
//    number of kmers in each block
//
//  Used as argument to merylFileReader::loadFromFile().
//  Populated by the file writer.

class merylFileIndex {
public:
  merylFileIndex() {
    clear();
  };

  void       set(kmpref  prefix,
                 FILE   *F,
                 uint64  nKmers) {

    if (_blockPosition == UINT64_MAX) {
      _blockPrefix   = prefix;
      _blockPosition = AS_UTL_ftell(F);
      _numKmers      = nKmers;
    }

    else {
      _numKmers     += nKmers;
    }

    if (_blockPrefix != prefix)
      fprintf(stderr, "set prefix 0x%s from [0x%s %8lu %8lu] to [%8lu %8lu]\n",
              toHex(prefix),
              toHex(_blockPrefix), _blockPosition, _numKmers,
              AS_UTL_ftell(F), _numKmers + nKmers);
    assert(_blockPrefix   == prefix);

    assert(_blockPosition <= AS_UTL_ftell(F));
  };

  void       clear(void) {
    _blockPrefix   = 0;
    _blockPosition = UINT64_MAX;
    _numKmers      = 0;
  }

  kmpref     blockPrefix(void)     { return((kmpref)_blockPrefix);   };
  uint64     blockPosition(void)   { return(        _blockPosition); };
  uint64     numKmers(void)        { return(        _numKmers);      };

private:
  uint64    _blockPrefix;     //  For compatibility, and alignment, _blockPrefix
  uint64    _blockPosition;   //  needs to be uint64 instead of the more correct
  uint64    _numKmers;        //  kmpref.
};



#endif  //  MERYL_UTIL_KMER_FILES_H
