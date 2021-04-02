
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

#ifndef MERYL_H
#define MERYL_H

#include "runtime.H"

#include "merylInput.H"
#include "merylOp.H"
#include "merylCountArray.H"

#include <stack>
#include <vector>
#include <algorithm>


class merylCommandBuilder {
public:
  merylCommandBuilder();
  ~merylCommandBuilder();

  void    initialize(char *opt);
  bool    processOptions(void);

  void    terminateOperation(void);
  bool    processOperation(void);

  bool    isOutput(void);
  bool    isPrinter(void);

  bool    isMerylInput(void);
  bool    isCanuInput(std::vector<char *> &err);
  bool    isSequenceInput(void);

  void    finalize(void);
  void    spawnThreads(void);

  void    printTree(merylOperation *op, uint32 indent);

public:
  uint32           numOperations(void)           { return(_opList.size()); };
  merylOperation  *getOperation(uint32 i)        { return(_opList[i]);     };

  uint32           numRoots(void)                { return(_opRoot.size()); };
  merylOperation  *getRoot(uint32 r)             { return(_opList   [ _opRoot[r] ]);  };
  merylOperation  *getRoot(uint32 r, uint32 t)   { return(_thList[t][ _opRoot[r] ]);  };

private:
  uint32    _terminating = 0;

  uint32    _optStringLen = 0;
  char      _optString[FILENAME_MAX + 1] = {0};

  char      _inoutName[FILENAME_MAX + 1] = {0};   //  Generic input or output name.
  char      _indexName[FILENAME_MAX + 1] = {0};   //  'merylIndex' in a meryl DB.
  char      _sqInfName[FILENAME_MAX + 1] = {0};   //  'info' in a Canu seqStore.
  char      _sqRdsName[FILENAME_MAX + 1] = {0};   //  'reads' in a Canu seqStore.

  //  Input from merylDB.

  //  Input from Canu seqStore.
  uint32    _segment       = 1;
  uint32    _segmentMax    = 1;

  //  Input from FASTA or FASTQ file.
  bool      _doCompression     = false;

  //  Output to merylDB.
  bool      _isOutput    = false;

  //  Output to dump file.
  bool      _printACGTorder = false;
  bool      _isPrint        = false;



  uint64    _allowedMemory;    //  These are set in the constructor,
  uint32    _allowedThreads;   //  based on what we know at run time.

  //  _opStack is a stack of operations, used when constructing the tree of operations.
  //
  //  _opList is a list of operations.

  std::stack <merylOperation *>   _opStack;
  std::vector<merylOperation *>   _opList;
  merylOperation                **_thList[64] = { nullptr };   //  Mirrors opList

  std::vector<uint32>             _opRoot;
};


#endif  //  MERYL_H
