
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

#ifndef MERYL_UTIL_KMER_WRITER_BLOCK_H
#define MERYL_UTIL_KMER_WRITER_BLOCK_H

#ifndef MERYL_UTIL_KMER_H
#error "include kmers.H, not this."
#endif

//  Responsible for taking a block of sorted kmers
//  and writing them to the appropriate file.

class merylFileWriter;

class merylBlockWriter {
public:
  merylBlockWriter(merylFileWriter *writer);
  ~merylBlockWriter();

public:
  void    addBlock(kmpref prefix, uint64 nKmers, kmdata *suffixes, kmvalu *values);

  void    finishBatch(void);
  void    finish(void);

private:
  void    closeFileDumpIndex(uint32 oi, uint32 iteration=UINT32_MAX);
  void    mergeBatches(uint32 oi);

private:
  merylFileWriter       *_writer;
  char                   _outName[FILENAME_MAX+1];

  //  Encoding data (copied from the input merylFileWriter).

  uint32                 _prefixSize;

  uint32                 _suffixSize;
  kmdata                 _suffixMask;

  uint32                 _numFilesBits;
  uint32                 _numBlocksBits;
  uint64                 _numFiles;
  uint64                 _numBlocks;

  //  File data.

  FILE                 **_datFiles;
  merylFileIndex       **_datFileIndex;

  //  Kmer data and et cetera.

  uint32                 _iteration;
};

#endif  //   MERYL_UTIL_KMER_WRITER_BLOCK_H
