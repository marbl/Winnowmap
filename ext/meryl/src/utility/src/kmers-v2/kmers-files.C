
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

#include "kmers.H"

namespace merylutil::inline kmers::v2 {

merylFileBlockReader::merylFileBlockReader() {
  _data        = NULL;

  _blockPrefix = 0;
  _nKmers      = 0;
  _nKmersMax   = 0;

  _kCode       = 0;
  _unaryBits   = 0;
  _binaryBits  = 0;
  _k1          = 0;

  _cCode       = 0;
  _c1          = 0;
  _c2          = 0;

  _suffixes    = NULL;
  _values      = NULL;
  _labels      = NULL;
}


merylFileBlockReader::~merylFileBlockReader() {
  delete    _data;
  delete [] _suffixes;
  delete [] _values;
  delete [] _labels;
}


bool
merylFileBlockReader::loadKmerFileBlock(FILE *inFile, uint32 activeFile, uint32 activeIteration) {

  //  If _data exists, we've already loaded the block, but haven't used it yet.

  if (_data)
    return(true);

  //  Otherwise, allocate _data, read the block from disk.  If nothing loaded,
  //  return false.

  _data = new stuffedBits(inFile);

  _blockPrefix = 0;
  _nKmers      = 0;

  if (_data->getLength() == 0) {
    delete _data;
    _data = NULL;

    return(false);
  }

  //  Decode the header of _data, but don't process the kmers yet.

  uint64 m1    = _data->getBinary(64);
  uint64 m2    = _data->getBinary(64);

  //  Version merylDataFile00 is the original format.
  //  It was replaced on 8 Oct 2020 by merylDataFile01.
  if      ((m1 == 0x7461446c7972656dllu) &&   //  merylDat
           (m2 == 0x0a3030656c694661llu)) {   //  aFile00\n
    _blockPrefix = _data->getBinary(64);
    _nKmers      = _data->getBinary(64);

    _kCode       = _data->getBinary(8);
    _unaryBits   = _data->getBinary(32);
    _binaryBits  = _data->getBinary(32);
    _k1          = _data->getBinary(64);

    _cCode       = _data->getBinary(8);
    _c1          = _data->getBinary(64);
    _c2          = _data->getBinary(64);

    _lCode       = 0;
    _labelBits   = 0;
    _l1          = 0;
    _l2          = 0;
  }

  //  Version merylDataFile01 added support for labels.
  else if ((m1 == 0x7461446c7972656dllu) &&   //  merylDat
           (m2 == 0x0a3130656c694661llu)) {   //  aFile01\n
    _blockPrefix = _data->getBinary(64);
    _nKmers      = _data->getBinary(64);

    _kCode       = _data->getBinary(8);
    _unaryBits   = _data->getBinary(32);
    _binaryBits  = _data->getBinary(32);
    _k1          = _data->getBinary(64);

    _cCode       = _data->getBinary(8);
    _c1          = _data->getBinary(64);
    _c2          = _data->getBinary(64);

    _lCode       = _data->getBinary(8);
    _labelBits   = _data->getBinary(6);
    _l1          = _data->getBinary(58);
    _l2          = _data->getBinary(64);
  }

  else if (((m1 & 0xffffffffffffffffllu) == 0x7461446c7972656dllu) &&   //  Match all but the
           ((m2 & 0xfff0f0ffffffffffllu) == 0x0a3030656c694661llu)) {   //  the version numbers.
    char v1 = (m2 >> 40) & 0xff;
    char v2 = (m2 >> 48) & 0xff;

    fprintf(stderr, "This version of meryl supports only data file versions 00 and 01.\n");
    fprintf(stderr, "This is a version %c%c data file.\n", v1, v2);
    fprintf(stderr, "  m1 = 0x%016" F_X64P "\n", m1);
    fprintf(stderr, "  m2 = 0x%016" F_X64P "\n", m2);
    exit(1);
  }

  else {
    fprintf(stderr, "merylFileReader::nextMer()-- Magic number mismatch in activeFile " F_U32 " activeIteration " F_U32 ".\n", activeFile, activeIteration);
    fprintf(stderr, "merylFileReader::nextMer()-- Expected 0x7461446c7972656d got 0x%016" F_X64P "\n", m1);
    fprintf(stderr, "merylFileReader::nextMer()-- Expected 0x0a3.30656c694661 got 0x%016" F_X64P "\n", m2);
    exit(1);
  }


#ifdef SHOW_LOAD
  fprintf(stderr, "loadKmerFileBlock()-- file %u iter %u:\n", activeFile, activeIteration);
  fprintf(stderr, "    prefix     0x%016lx\n", _blockPrefix);
  fprintf(stderr, "    nKmers     " F_U64 "\n", _nKmers);
  fprintf(stderr, "    kCode      " F_U32 "\n", _kCode);
  fprintf(stderr, "    unaryBits  " F_U32 "\n", _unaryBits);
  fprintf(stderr, "    binaryBits " F_U32 "\n", _binaryBits);
  fprintf(stderr, "    k1         " F_U64 "\n", _k1);
  fprintf(stderr, "    cCode      " F_U32 "\n", _cCode);
  fprintf(stderr, "    c1         " F_U64 "\n", _c1);
  fprintf(stderr, "    c2         " F_U64 "\n", _c2);
  fprintf(stderr, "    lCode      " F_U32 "\n", _lCode);
  fprintf(stderr, "    labelBits  " F_U32 "\n", _labelBits);
  fprintf(stderr, "    l1         " F_U64 "\n", _l1);
  fprintf(stderr, "    l2         " F_U64 "\n", _l2);
#endif

  return(true);
}


#if 0
void
merylFileBlockReader::decodeBlock(void) {
  if (_data == NULL)
    return;

  resizeArrayPair(_suffixes, _values, 0, _nKmersMax, _nKmers, _raAct::doNothing);
  decodeBlock(_suffixes, _values);
}
#endif


void
merylFileBlockReader::decodeKmerFileBlockData(kmdata *suffixes) {
  if      (_kCode == 1) {
    kmdata  thisPrefix = 0;

    for (uint32 kk=0; kk<_nKmers; kk++) {
      thisPrefix += (kmdata)_data->getUnary();

      uint32 ls = (_binaryBits <= 64) ? (0)           : (_binaryBits - 64);
      uint32 rs = (_binaryBits <= 64) ? (_binaryBits) : (64);

      suffixes[kk]   = thisPrefix;
      suffixes[kk] <<= ls;
      suffixes[kk]  |= _data->getBinary(ls);
      suffixes[kk] <<= rs;
      suffixes[kk]  |= _data->getBinary(rs);
    }
  }

  else {
    fprintf(stderr, "ERROR: unknown kCode 0x%02x\n", _kCode), exit(1);
  }
}


void
merylFileBlockReader::decodeKmerFileBlockValu(kmvalu *values) {
  if      (_cCode == 0) {
    for (uint32 kk=0; kk<_nKmers; kk++)
      values[kk] = 0;
  }

  else if (_cCode == 1) {
    for (uint32 kk=0; kk<_nKmers; kk++)
      values[kk] = _data->getBinary(32);
  }

  else if (_cCode == 2) {
    for (uint32 kk=0; kk<_nKmers; kk++)
      values[kk] = _data->getBinary(64);
  }

  else {
    fprintf(stderr, "ERROR: unknown cCode 0x%02x\n", _cCode), exit(1);
  }
}


void
merylFileBlockReader::decodeKmerFileBlockLabl(kmlabl *labels) {
  if      (_lCode == 0) {
    for (uint32 kk=0; kk<_nKmers; kk++)
      labels[kk] = 0;
  }

  else if (_lCode == 1) {
    for (uint32 kk=0; kk<_nKmers; kk++)
      labels[kk] = _data->getBinary(_labelBits);
  }

  else {
    fprintf(stderr, "ERROR: unknown lCode 0x%02x\n", _lCode), exit(1);
  }
}



//  Decode a block of kmers.  The block has _nKmers kmers in it, but the
//  arrays have allocated only _nKmersMax space.
//
void
merylFileBlockReader::decodeKmerFileBlock(void) {

  if (_data == nullptr)
    return;

  resizeArray(_suffixes, _values, _labels, 0, _nKmersMax, _nKmers, _raAct::doNothing);

  decodeKmerFileBlockData(_suffixes);
  decodeKmerFileBlockValu(_values);
  decodeKmerFileBlockLabl(_labels);

  delete _data;   _data = nullptr;
}


//  If there is data to decode, decode it (if space supplied to decode into)
//  the delete the raw data.
//
void
merylFileBlockReader::decodeKmerFileBlock(kmdata *suffixes, kmvalu *values, kmlabl *labels) {

  if (_data == nullptr)
    return;

  if (suffixes)   decodeKmerFileBlockData(suffixes);
  if (values)     decodeKmerFileBlockValu(values);
  if (labels)     decodeKmerFileBlockLabl(labels);

  delete _data;   _data = nullptr;
}

}  //  namespace merylutil::kmers::v2

