
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

#ifndef MERYLOPCOMPUTE_H
#define MERYLOPCOMPUTE_H

#ifndef MERYLINCLUDE
#error "Do not use merylOpCompute.H, use meryl.H instead."
#endif


//  This is used to hold the input file index, the value and the label from
//  each input.  Two lists are created, one of the inputs with the smallest
//  kmer (_act) and one for all inputs (_inp).
//
//  For _act, _idx is the index of the input.
//  For _inp, _idx is 0 if this is input has the smallest kmer,
//                    uint32max otherwise.
//
struct merylActList {
  uint32  _idx;   //  Index of the input for _act; 0 if kmer is valid in _inp
  kmvalu  _val;   //  Value of this kmer
  kmlabl  _lab;   //  Label of this kmer
};

constexpr
uint32
merylActListMax = 128;



//  This class is used for performing the compute.
//
//  It is called in merylCommandBuilder::spawnThreads().  The
//  compute is done using these objects.
//
class merylOpCompute {
public:
  merylOpCompute(merylOpTemplate *op, uint32 slice, uint32 nInputs);
  ~merylOpCompute();

  //  Methods used to run this operation.
public:
  bool                  nextMer(void);
  kmer                  theFMer(void)            { return(_kmer);   };

private:
  void                  findOutputKmer(void);
  void                  findOutputValue(void);
  void                  findOutputLabel(void);

  //
  //  Methods to set up this operations.
  //
public:
  void                  addSliceInput(merylInput *i) {
    _inputs.push_back(i);
  }

  void                  addSliceOutput    (merylFileWriter *outDbse, uint32 slice);
  void                  addSliceLister    (opPname which, char const *listName, compressedFileWriter *list, uint32 slice);
  void                  addSliceStatistics(merylOpTemplate *ot, uint32 slice);

  //
  //  Pointer back to the template, so we can get other useful bits instead
  //  of copying everything here.
  //

  merylOpTemplate               *_ot   = nullptr;

  //
  //  Inputs.  A list of the inputs to this operation, and the specific slice
  //  we're working on.  The input objects are all new, different from those
  //  in the template, so we can limit each one to a specific database slice.
  //

  std::vector<merylInput *>      _inputs;

  //
  //  Per-slice ANALYSES.
  //

  merylHistogram                *_statsAcc       = nullptr;    //  Per-thread stats

  //
  //  Per-slice OUTPUTS.
  //

  merylStreamWriter             *_outDbseSlice   = nullptr;    //  Per-thread output, I own this.
  merylStreamWriter             *_outDbse        = nullptr;    //  Always _outDbseSlice

  compressedFileWriter          *_outListSlice   = nullptr;    //  Per-thread output, I own this.
  compressedFileWriter          *_outList        = nullptr;    //  Either _outListSlice or ot->_outList

  //mpressedFileWriter          *_outShowSlice   = nullptr;    //  Unused.
  compressedFileWriter          *_outShow        = nullptr;    //  Always ot->_outShow

  //ar const                    *_outPipe        = nullptr;

  //
  //  Data used to compute the action.
  //

  kmer                           _kmer;               //  The kmer we're computing.

  uint32                         _actLen = 0;

  merylActList                   _acts[merylActListMax];
  merylActList                   _inps[merylActListMax];

  merylActList                  *_actAlloc;
  merylActList                  *_acta  = nullptr;
  merylActList                  *_inpa  = nullptr;

  //
  //  Selection.
  //

  bool   shouldKmerBeOutput(uint32 ii);
  bool   shouldKmerBeOutput(void);

  friend class merylCommandBuilder;
};


#endif  //  MERYLOPCOMPUTE_H
