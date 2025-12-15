
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

#ifndef MERYLUTIL_DS_ARRAYS_V1_H
#define MERYLUTIL_DS_ARRAYS_V1_H

#include "types.H"

namespace merylutil::inline arrays::inline v1 {

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



//  Return a copy of the NUL terminated string 'fr'.
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



//  Copy NUL terminated string 'fr' into 'to' (re)allocating
//  space as needed.  Return the 'to' pointer.
template<typename TT, typename LL>
TT *
duplicateString(TT const *fr, TT *&to, LL &toLen, LL &toMax) {

  if (fr == nullptr) {
    toLen = 0;
    increaseArray(to, toLen+1, toMax, 16, _raAct::doNothing);
    to[0] = 0;
  }
  else {
    toLen = strlen(fr);
    increaseArray(to, toLen+1, toMax, 16, _raAct::doNothing);
    memcpy(to, fr, sizeof(TT) * (toLen+1));
  }

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



//
//  Change the size of an array.  The array is ALWAYS reallocated.
//
//  Elements        0 .. arrayLen are copied from the old to the new array.
//  Elements arrayLen .. newMax   can be (optionally) set to zero.
//
//  WARNING: If arrayLen is larger than newMax, only the first newMax
//           elements are copied.  arrayLen IS NOT modified to reflect this.
//
//  Output parameter arrayMax is set to newMax; arrayMax is not otherwise used.
//
//  WARNING: Use increaseArray() or resizeArray() instead of this function.
//
template<typename TT, typename LL>
void
setArraySize(TT*& array, uint64 arrayLen, LL &arrayMax, uint64 newMax, _raAct op=_raAct::copyData, bool update=true) {

  arrayLen = std::min(newMax, arrayLen);

  TT *copy = new TT [newMax];

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
    for (uint32 ii=arrayLen; ii<newMax; ii++)
      copy[ii] = TT();

  if (update)
    arrayMax = newMax;
}



//
//  Ensure that there is enough space to hold one more element in the array.
//  Increase the array by 'moreSpace' elements if needed.
//
//  With the array used as a stack, where elements are added to the tail:
//
//    while (...) {
//      increaseArray(arr, arrLen, arrMax, 32)
//      arr[arrLen++] = ...
//    }
//
//  will allocate 32 more elements to the array when arrLen == arrMax, and do
//  nothing otherwise.  After the call, array element arr[arrLen] is
//  guaranteed to exist.
//
//  Likewise, for random-access, element arr[idx] is guaranteed to exist
//  after the call:
//
//    idx = arrMax + ...
//    increaseArray(arr, idx, arrMax, 32)
//    arr[idx] = ...
//
//  will do nothing if idx < arrMax, but will resize the array to have idx+32
//  elements otherwise.  In this case, the 'length' of the array is not
//  tracked, we only know that the array has arrMax allocated elements (and
//  you probably want to set op=_raAct::copyDataClearNew to initialize any
//  newly allocated elements to zero).
//  
//  If the array is reallocated, the contents of the entire array (0
//  .. arrayMax) are copied to the new space.  By default, new elements are
//  NOT cleared to zero; set 'op' to '_raAct::copyDataClearNew' to enable.
//
//  'moreSpace' must be at least 1.
//

static   //  Helper for increaseArray() -- return the first muliple of 'moreSpace'
inline   //  that is strictly larger than 'idx'.
uint64
nmc(uint64 idx, uint64 moreSpace) {
  uint64 ns = (moreSpace <= 1) ? (            (idx            ) + 1)
                               : (moreSpace * (idx / moreSpace) + moreSpace);
  assert(idx < ns);

  return ns;
}

template<typename TT, typename LL>
void
increaseArray(TT*& array, uint64 idx, LL &arrayMax, uint64 moreSpace, _raAct op=_raAct::copyData) {
  if (idx < arrayMax)
    return;
  setArraySize(array, arrayMax, arrayMax, nmc(idx, moreSpace), op);
}


template<typename T1, typename T2, typename LL>
void
increaseArrayPair(T1*& array1, T2*& array2, uint64 idx, LL &arrayMax, uint64 moreSpace, _raAct op=_raAct::copyData) {
  if (idx < arrayMax)
    return;
  setArraySize(array1, arrayMax, arrayMax, nmc(idx, moreSpace), op, false);
  setArraySize(array2, arrayMax, arrayMax, nmc(idx, moreSpace), op);
}


template<typename T1, typename T2, typename LL>
void
increaseArray(T1*& array1, T2*& array2, uint64 idx, LL &arrayMax, uint64 moreSpace, _raAct op=_raAct::copyData) {
  if (idx < arrayMax)
    return;
  setArraySize(array1, arrayMax, arrayMax, nmc(idx, moreSpace), op, false);
  setArraySize(array2, arrayMax, arrayMax, nmc(idx, moreSpace), op);
}


template<typename T1, typename T2, typename T3, typename LL>
void
increaseArray(T1*& array1, T2*& array2, T3*& array3, uint64 idx, LL &arrayMax, uint64 moreSpace, _raAct op=_raAct::copyData) {
  if (idx < arrayMax)
    return;
  setArraySize(array1, arrayMax, arrayMax, nmc(idx, moreSpace), op, false);
  setArraySize(array2, arrayMax, arrayMax, nmc(idx, moreSpace), op, false);
  setArraySize(array3, arrayMax, arrayMax, nmc(idx, moreSpace), op);
}


template<typename T1, typename T2, typename T3, typename T4, typename LL>
void
increaseArray(T1*& array1, T2*& array2, T3*& array3, T3*& array4, uint64 idx, LL &arrayMax, uint64 moreSpace, _raAct op=_raAct::copyData) {
  if (idx < arrayMax)
    return;
  setArraySize(array1, arrayMax, arrayMax, nmc(idx, moreSpace), op, false);
  setArraySize(array2, arrayMax, arrayMax, nmc(idx, moreSpace), op, false);
  setArraySize(array3, arrayMax, arrayMax, nmc(idx, moreSpace), op, false);
  setArraySize(array4, arrayMax, arrayMax, nmc(idx, moreSpace), op);
}



//
//  Ensure the array has at least 'newMax' elements:
//    if arrayMax < newMax  ->  resize to newMax exactly.
//    else                  ->  do nothing.
//
//  By default, new elements are not initialized, use
//  _raAct::copyDataClearNew if desired.
//
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


template<typename T1, typename T2, typename LL>
void
resizeArray(T1*& array1, T2*& array2, uint64 arrayLen, LL &arrayMax, uint64 newMax, _raAct op=_raAct::copyData) {

  if (newMax <= arrayMax)
    return;

  setArraySize(array1, arrayLen, arrayMax, newMax, op);
  setArraySize(array2, arrayLen, arrayMax, newMax, op);
}


template<typename T1, typename T2, typename T3, typename LL>
void
resizeArray(T1*& array1, T2*& array2, T3*& array3, uint64 arrayLen, LL &arrayMax, uint64 newMax, _raAct op=_raAct::copyData) {

  if (newMax <= arrayMax)
    return;

  setArraySize(array1, arrayLen, arrayMax, newMax, op);
  setArraySize(array2, arrayLen, arrayMax, newMax, op);
  setArraySize(array3, arrayLen, arrayMax, newMax, op);
}


template<typename T1, typename T2, typename T3, typename T4, typename LL>
void
resizeArray(T1*& array1, T2*& array2, T3*& array3, T4*& array4, uint64 arrayLen, LL &arrayMax, uint64 newMax, _raAct op=_raAct::copyData) {

  if (newMax <= arrayMax)
    return;

  setArraySize(array1, arrayLen, arrayMax, newMax, op);
  setArraySize(array2, arrayLen, arrayMax, newMax, op);
  setArraySize(array3, arrayLen, arrayMax, newMax, op);
  setArraySize(array4, arrayLen, arrayMax, newMax, op);
}



////////////////////////////////////////
//
//  Append an element to an array.
//
template<typename T, typename L>
L
appendToArray(T elt, T*& array, L &arrayLen, L &arrayMax, uint64 moreSpace=512, _raAct op=_raAct::copyData) {
  increaseArray(array, arrayLen, arrayMax, moreSpace, op);
  array[arrayLen++] = elt;
  return arrayLen - 1;
}


}  //  merylutil::arrays::v1

#endif  //  MERYLUTIL_DS_ARRAYS_V1_H
