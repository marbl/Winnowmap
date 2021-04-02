
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

#ifndef MERYL_UTIL_KMER_WRITER_STREAM_H
#define MERYL_UTIL_KMER_WRITER_STREAM_H

#ifndef MERYL_UTIL_KMER_H
#error "include kmers.H, not this."
#endif


class merylFileWriter;


class merylStreamWriter {
public:
  merylStreamWriter(merylFileWriter *writer, uint32 fileNumber);
  ~merylStreamWriter();

public:
  void    addMer(kmer k, kmvalu c);

private:
  void    dumpBlock(kmpref nextPrefix=~((kmpref)0));

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

  uint32                 _filePrefix;

  FILE                  *_datFile;
  merylFileIndex        *_datFileIndex;

  //  Kmer data and etc for writing the stream.

  kmpref                 _batchPrefix;
  uint64                 _batchNumKmers;
  uint64                 _batchMaxKmers;
  kmdata                *_batchSuffixes;
  kmvalu                *_batchValues;
};

#endif  //   MERYL_UTIL_KMER_WRITER_STREAM_H
