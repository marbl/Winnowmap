
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

#ifndef MERYL_UTIL_KMER_FILES_WRITER_H
#define MERYL_UTIL_KMER_FILES_WRITER_H

#ifndef MERYL_UTIL_KMER_H
#error "include kmers.H, not this."
#endif

#include "kmers-files.H"
#include "kmers-writer-block.H"
#include "kmers-writer-stream.H"

//  Responsible for writing data to a file.
//
//  It's a little confusing to use; the merylFileWriter is used to get a
//  blockWriter or streamWriter object that actually accepts the data to
//  write.

class merylFileWriter {
public:
  merylFileWriter(const char *outputName,
                      uint32      prefixSize = 0);

  ~merylFileWriter();

public:
  char   *filename(void)  { return(_outName);  };

public:
  void    initialize(uint32 prefixSize = 0, bool isMultiSet = false);

  merylBlockWriter  *getBlockWriter(void)        { return(new merylBlockWriter (this));      };
  merylStreamWriter *getStreamWriter(uint32 ff)  { return(new merylStreamWriter(this, ff));  };

public:
  uint32  numberOfFiles(void)           { return(_numFiles);                      };
  uint64  firstPrefixInFile(uint32 ff)  { return(((uint64)ff) << _numBlocksBits); };
  uint64  lastPrefixInFile(uint32 ff)   { return(firstPrefixInFile(ff + 1) - 1);  };

  uint32  fileNumber(uint64 prefix);

private:
  void    writeBlockToFile(FILE            *datFile,
                           merylFileIndex  *datFileIndex,
                           kmpref           blockPrefix,
                           uint64           nKmers,
                           kmdata          *suffixes,
                           kmvalu          *values);

private:
  bool                       _initialized;

  char                       _outName[FILENAME_MAX+1];

  uint32                     _prefixSize;

  uint32                     _suffixSize;
  kmdata                     _suffixMask;

  uint32                     _numFilesBits;
  uint32                     _numBlocksBits;
  uint64                     _numFiles;
  uint64                     _numBlocks;

  bool                       _isMultiSet;

  merylHistogram             _stats;

  friend class merylBlockWriter;
  friend class merylStreamWriter;
};




#endif  //  MERYL_UTIL_KMER_FILES_WRITER_H
