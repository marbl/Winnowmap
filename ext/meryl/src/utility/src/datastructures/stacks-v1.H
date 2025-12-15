
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-JUL-20
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef MERYLUTIL_DATASTRUCTURES_STACKS_V1_H
#define MERYLUTIL_DATASTRUCTURES_STACKS_V1_H

#include "types.H"

namespace merylutil::inline stacks::inline v1 {


template <typename T>
class stack {
public:
  stack()   {               }
  ~stack()  { delete [] _S; }

  uint64    depth(void)     { return _Stop;      }
  bool      empty(void)     { return _Stop == 0; }

  stack<T> &clear(void)     { _Stop = 0;  return *this; }

  T         pop(void);
  T        &top(void);

  T        &operator[](uint64 i);

  T         push(T const item);
  stack<T> &exchange(void);


private:
  uint64  _Smax = 0;
  uint64  _Stop = 0;
  T      *_S    = nullptr;
};


template <typename T>
T
stack<T>::pop(void) {
  if (_Stop < 1)
    fprintf(stderr, "ERROR: stack::pop() has no elements to return.\n");
  assert(_Stop >= 1);
  return _S[--_Stop];
}


template <typename T>
T&
stack<T>::top(void) {
  if (_Stop < 1)
    fprintf(stderr, "ERROR: stack::top() has no elements to return.\n");
  assert(_Stop >= 1);
  return _S[_Stop-1];
}

template <typename T>
T&
stack<T>::operator[](uint64 i) {
  if (_Stop <= i)
    fprintf(stderr, "ERROR: stack::operator[]() has only %lu elements, cannot access element %lu.\n", _Stop, i);
  assert(i < _Stop);
  return _S[i];
}


template <typename T>
T
stack<T>::push(T const item) {
  increaseArray(_S, _Stop, _Smax, 128);
  return _S[_Stop++] = item;
}


template <typename T>
stack<T> &
stack<T>::exchange(void) {
  if (_Stop < 2)
    fprintf(stderr, "ERROR: stack::exchange() doesn't have enough elements to exchange.\n");
  assert(_Stop >= 2);
  std::swap(_S[_Stop-2], _S[_Stop-1]);
  return *this;
}


}  //    merylutil::stacks::v1

#endif  //  MERYLUTIL_DATASTRUCTURES_STACKS_V1_H
