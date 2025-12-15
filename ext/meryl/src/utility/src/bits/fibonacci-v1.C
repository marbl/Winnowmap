
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

#include "types.H"
#include "fibonacci-v1.H"

//  Set up the Fibonacci encoding table.
//
//  It takes 46 values to saturate a uint32 (fib[47] > uint32).
//  It takes 92 values to saturate a uint64 (fib[93] > uint64).
//
//  This is NOT an efficient method; it should be used only to get a
//  compile-time evaluatable function for populating fibonaciNumber[].

namespace merylutil::inline bits::inline v1 {

constexpr
uint64
__fibNum(uint32 f) {
  uint64  fs[93] = { 1, 1, 0 };

  for (uint32 ii=2; ii<=f; ii++)
    fs[ii] = fs[ii-1] + fs[ii-2];

  return(fs[f]);
}  //  C++14 required!

uint64
__fibNumV[93] = {
  __fibNum( 0), __fibNum( 1), __fibNum( 2), __fibNum( 3), __fibNum( 4), __fibNum( 5), __fibNum( 6), __fibNum( 7), __fibNum( 8), __fibNum( 9),
  __fibNum(10), __fibNum(11), __fibNum(12), __fibNum(13), __fibNum(14), __fibNum(15), __fibNum(16), __fibNum(17), __fibNum(18), __fibNum(19),
  __fibNum(20), __fibNum(21), __fibNum(22), __fibNum(23), __fibNum(24), __fibNum(25), __fibNum(26), __fibNum(27), __fibNum(28), __fibNum(29),
  __fibNum(30), __fibNum(31), __fibNum(32), __fibNum(33), __fibNum(34), __fibNum(35), __fibNum(36), __fibNum(37), __fibNum(38), __fibNum(39),
  __fibNum(40), __fibNum(41), __fibNum(42), __fibNum(43), __fibNum(44), __fibNum(45), __fibNum(46), __fibNum(47), __fibNum(48), __fibNum(49),
  __fibNum(50), __fibNum(51), __fibNum(52), __fibNum(53), __fibNum(54), __fibNum(55), __fibNum(56), __fibNum(57), __fibNum(58), __fibNum(59),
  __fibNum(60), __fibNum(61), __fibNum(62), __fibNum(63), __fibNum(64), __fibNum(65), __fibNum(66), __fibNum(67), __fibNum(68), __fibNum(69),
  __fibNum(70), __fibNum(71), __fibNum(72), __fibNum(73), __fibNum(74), __fibNum(75), __fibNum(76), __fibNum(77), __fibNum(78), __fibNum(79),
  __fibNum(80), __fibNum(81), __fibNum(82), __fibNum(83), __fibNum(84), __fibNum(85), __fibNum(86), __fibNum(87), __fibNum(88), __fibNum(89),
  __fibNum(90), __fibNum(91), __fibNum(92)
};

static_assert(__fibNum(45) < __fibNum(46));  //  Fails if 32-bit signed
static_assert(__fibNum(46) < __fibNum(47));  //  Fails if 32-bit unsigned
static_assert(__fibNum(91) < __fibNum(92));  //  Fails if 64-bit signed
//atic_assert(__fibNum(92) < __fibNum(93));  //  Always fails.

}  //  namespace merylutil::bits::v1
