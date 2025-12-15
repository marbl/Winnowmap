
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

#ifndef MERYLUTIL_KMERS_READER_V2_H
#define MERYLUTIL_KMERS_READER_V2_H

#ifndef MERYLUTIL_KMERS_H
#error "include kmers.H, not this."
#endif

namespace merylutil::inline kmers::v2 {

class merylFileReader {
private:
  stuffedBits  *openMasterIndex(void);

  void          initializeFromMasterI_v00(void);
  void          initializeFromMasterI_v01(stuffedBits  *masterIndex);
  void          initializeFromMasterI_v02(stuffedBits  *masterIndex);
  void          initializeFromMasterI_v03(stuffedBits  *masterIndex);
  void          initializeFromMasterI_v04(stuffedBits  *masterIndex);
  bool          initializeFromMasterIndex(std::vector<char const *> *errors = nullptr);

public:
  merylFileReader(const char *inputName, std::vector<char const *> *errors=nullptr);

  merylFileReader(const char *inputName,
                  uint32      threadFile, std::vector<char const *> *errors=nullptr);

  ~merylFileReader();

public:
  void    rewind(void) {
    _activeMer  = 0;    //  Position we are at in the block loaded.
    _activeFile = 0;

    _nKmers     = 0;    //  Kmers in the block we have loaded.

    if (_threadFile != UINT32_MAX)
      _activeFile = _threadFile;

    merylutil::closeFile(_datFile);
  };

public:
  void    loadStatistics(std::vector<char const *> *errors = nullptr);
  void    dropStatistics(void);

public:
  void    enableThreads(uint32 threadFile);

public:
  void    loadBlockIndex(void);

public:
  bool    nextMer(void);

  kmer    theFMer(void)        { return(_kmer);        };

  //#warning OBSOLETE theValue() and theLabel()
  kmvalu  theValue(void)       { return(_kmer._val);   };
  kmlabl  theLabel(void)       { return(_kmer._lab);   };

  //bool    isMultiSet(void)     { return(_isMultiSet);  };

  char   *filename(void)       { return(_inName);      };

  merylHistogram  *stats(void) {
    loadStatistics();
    return(_stats);
  }

  //  For direct access to the kmer blocks.
public:
  uint32  prefixSize(void)     { return(_prefixSize); };
  uint32  suffixSize(void)     { return(_suffixSize); };

  uint32  numFiles(void)       { return(_numFiles);  };
  uint32  numBlocks(void)      { return(_numBlocks); };

  uint32  numFilesBits(void)   { return(_numFilesBits);  };
  uint32  numBlocksBits(void)  { return(_numBlocksBits); };

  FILE   *blockFile(uint32 ff) {
    FILE *F = NULL;

    if (ff < _numFiles)
      F = openInputBlock(_inName, ff, _numFiles);

    return(F);
  };

  merylFileIndex   &blockIndex(uint32 bb) {
    return(_blockIndex[bb]);
  };

private:
  char                       _inName[FILENAME_MAX+1] = {0};

  uint32                     _prefixSize    = 0;
  uint32                     _suffixSize    = 0;
  uint32                     _numFilesBits  = 0;
  uint32                     _numBlocksBits = 0;

  uint32                     _numFiles      = 0;
  uint32                     _numBlocks     = 0;

  bool                       _isMultiSet    = false;

  merylHistogram            *_stats         = nullptr;
  uint32                     _statsVersion  = 0;
  uint64                     _statsOffset   = 0;

  FILE                      *_datFile       = nullptr;

  merylFileBlockReader      *_block         = nullptr;
  merylFileIndex            *_blockIndex    = nullptr;

  kmer                       _kmer          = kmer();

  uint64                     _prefix        = 0;

  kmdata                     _activeMer     = 0;
  uint32                     _activeFile    = 0;

  uint32                     _threadFile    = UINT32_MAX;

  uint64                     _nKmers        = 0;
  uint64                     _nKmersMax     = 0;
  kmdata                    *_suffixes      = nullptr;
  kmvalu                    *_values        = nullptr;
  kmlabl                    *_labels        = nullptr;
};

}  //  namespace merylutil::kmers::v2

#endif  //  MERYLUTIL_KMERS_READER_V2_H
