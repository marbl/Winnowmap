
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

#ifndef ARRAYS_H
#define ARRAYS_H

#include "types.H"
#include <algorithm>


enum class _raAct {
  doNothing        = 0x00,
  copyData         = 0x01,
  clearNew         = 0x02,
  copyDataClearNew = 0x03,
};


inline   //  Combine two _raAct into one.
_raAct
operator|(_raAct a, _raAct b) {

  if (a == _raAct::doNothing)  return(b);
  if (b == _raAct::doNothing)  return(a);

  if ((a == _raAct::copyData) && (b == _raAct::copyData))   return(_raAct::copyData);
  if ((a == _raAct::copyData) && (b == _raAct::clearNew))   return(_raAct::copyDataClearNew);
  if ((a == _raAct::clearNew) && (b == _raAct::copyData))   return(_raAct::copyDataClearNew);
  if ((a == _raAct::clearNew) && (b == _raAct::clearNew))   return(_raAct::clearNew);

  if (a == _raAct::copyDataClearNew)  return(_raAct::copyDataClearNew);
  if (b == _raAct::copyDataClearNew)  return(_raAct::copyDataClearNew);

  assert(0);
  return(_raAct::doNothing);
}


inline   //  Return true if _raAct a has property b set.
bool
operator&(_raAct a, _raAct b) {

  if ((a == _raAct::copyData)         && (b == _raAct::copyData))           return(true);
  if ((a == _raAct::copyDataClearNew) && (b == _raAct::copyData))           return(true);

  if ((a == _raAct::clearNew)         && (b == _raAct::clearNew))           return(true);
  if ((a == _raAct::copyDataClearNew) && (b == _raAct::clearNew))           return(true);

  if ((a == _raAct::copyDataClearNew) && (b == _raAct::copyDataClearNew))   return(true);

  return(false);
}


//  Allocate an array of size 'allocSize', and set 'arrayMax' to that value.
//  By default. clear the array.
template<typename TT, typename LL>
void
allocateArray(TT*& array, LL &arrayMax, uint64 allocSize, _raAct op=_raAct::clearNew) {

  if (array != NULL)
    delete [] array;

  arrayMax = allocSize;
  array    = new TT [allocSize];

  assert(arrayMax == allocSize);   //  Make sure we don't truncate the value!

  if (op == _raAct::clearNew)
    memset(array, 0, sizeof(TT) * allocSize);
}


//  Allocate an array of size 'allocSize'.
//  By default, clear the array.
template<typename TT>
void
allocateArray(TT*& array, uint64 allocSize, _raAct op=_raAct::clearNew) {

  if (array != NULL)
    delete [] array;

  array    = new TT [allocSize];

  if (op == _raAct::clearNew)
    memset(array, 0, sizeof(TT) * allocSize);
}



template<typename TT>
TT *
duplicateString(TT const *fr) {

  if (fr == NULL)
    return(NULL);

  uint32  ln = strlen(fr);
  TT     *to = new TT [ln+1];

  memcpy(to, fr, sizeof(TT) * (ln+1));

  return(to);
}



template<typename TT, typename LL>
void
duplicateArray(TT*& to, LL &toLen, LL &toMax, TT const *fr, LL frLen, LL frMax=0, bool forceAlloc=false) {

  if (fr == NULL)
    assert(frLen == 0);

  if ((toMax < frLen) || (forceAlloc)) {
    delete [] to;

    toMax = frLen;
    to    = new TT [toMax];
  }

  toLen = frLen;

  if (frLen > 0)
    memcpy(to, fr, sizeof(TT) * frLen);
}


//  Set the array size to 'newMax'.
//  No guards, the array will ALWAYS be reallocated.
//
template<typename TT, typename LL>
void
setArraySize(TT*& array, uint64 arrayLen, LL &arrayMax, uint64 newMax, _raAct op=_raAct::copyData) {

  arrayMax =          newMax;
  arrayLen = std::min(newMax, arrayLen);

  TT *copy = new TT [arrayMax];

  if ((array != nullptr) &&
      (arrayLen > 0) &&
      ((op == _raAct::copyData) ||
       (op == _raAct::copyDataClearNew)))
    for (uint32 ii=0; ii<arrayLen; ii++)
      copy[ii] = array[ii];

  delete [] array;
  array = copy;

  if ((op == _raAct::clearNew) ||
      (op == _raAct::copyDataClearNew))
    for (uint32 ii=arrayLen; ii<arrayMax; ii++)
      copy[ii] = TT();
}



//  Ensure that there is enough space to hold one more element in the array.
//  Increase the array by 'moreSpace' if needed.
//
//  With the array used as a stack, a call of
//    increaseArray(arr, arrLen, arrMax, 32)
//  will allocate 32 more elements if arrLen == arrMax, and do nothing
//  otherwise.  After the call, array element arr[arrLen] is guaranteed to
//  exist.  If arrLen > arrMax, see below.
//
//  With the array used for random access, a call of
//    increaseArray(arr, idx, arrMax, 32)
//  will do nothing if idx < arrMax, and resize the array to have idx+32
//  elements otherwise.
//  
//  In both cases, if 'moreSpace' is 0, it is reset to 1.
//
//  If the array is reallocated, the contents of the entire array are copied
//  to the new space.  New elements are NOT cleared to zero; override op as
//  desired.
//
template<typename TT, typename LL>
void
increaseArray(TT*& array, uint64 idx, LL &arrayMax, uint64 moreSpace, _raAct op=_raAct::copyData) {
  uint64  newMax = idx + ((moreSpace == 0) ? 1 : moreSpace);

  if (idx < arrayMax)
    return;

  setArraySize(array, arrayMax, arrayMax, newMax, op);
}


template<typename T1, typename T2, typename LL>
void
increaseArrayPair(T1*& array1, T2*& array2, uint64 idx, LL &arrayMax, uint64 moreSpace, _raAct op=_raAct::copyData) {
  uint64  newMax = idx + ((moreSpace == 0) ? 1 : moreSpace);

  if (idx < arrayMax)
    return;

  setArraySize(array1, arrayMax, arrayMax, newMax, op);
  setArraySize(array2, arrayMax, arrayMax, newMax, op);
}


//  Resize the array so that it is at least as big as new max.  Do nothing
//  if the array is big enough already.

template<typename TT, typename LL>
void
resizeArray(TT*& array, uint64 arrayLen, LL &arrayMax, uint64 newMax, _raAct op=_raAct::copyData) {

  if (newMax <= arrayMax)
    return;

  setArraySize(array, arrayLen, arrayMax, newMax, op);
}


template<typename T1, typename T2, typename LL>
void
resizeArrayPair(T1*& array1, T2*& array2, uint64 arrayLen, LL &arrayMax, uint64 newMax, _raAct op=_raAct::copyData) {

  if (newMax <= arrayMax)
    return;

  setArraySize(array1, arrayLen, arrayMax, newMax, op);
  setArraySize(array2, arrayLen, arrayMax, newMax, op);
}


#endif  //  ARRAYS_H
