
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

#ifndef MERYLUTIL_SEQUENCE_V1_H
#define MERYLUTIL_SEQUENCE_V1_H

#include "types.H"
#include "files.H"

namespace merylutil::inline sequence::inline v1 {

//////////
//
//  Various functions for reverse complementing sequences.
//    reverseComplementSequence() is in place.
//    reverseComplementCopy()     returns a 'new' character array of length 'len+1'.
//    reverseComplement()         reverses both sequence and bases; used in tgTig still.
//
void   reverseComplementSequence(char *seq, int len);
char  *reverseComplementCopy    (char *seq, int len);

template<typename qvType>
void   reverseComplement(char *seq, qvType *qlt, int len);


//////////
//
//  Compress homopolymer runs to a single letter.  More comments in the .C file.
//
uint32 homopolyCompress(char   *bases, uint32 basesLen,
                        char   *compr=nullptr,
                        uint32 *ntoc=nullptr, char skip=0x7f);

//////////
//
//  Encode/decode sequence between ASCII and binary.
//    Encoding accepts upper or lower case letters.
//    Decoding ALWAYS returns uppercase letters (c=0).
//
//  The two and three bit encodings map ACGTN to 01324 and relies on the
//  ASCII encoding (where 'c' denotes uppercase (c=0) or lowercase (c=1)):
//      A      a     01c0 000 1 == '0' * 4 == 0  -> 0x...............0 -> 0
//      C      c     01c0 001 1 == '1' * 4 == 4  -> 0x..............1. -> 1
//      T      t     01c1 010 0 == '2' * 4 == 8  -> 0x.............2.. -> 2
//      G      g     01c0 011 1 == '3' * 4 == 12 -> 0x............3... -> 3
//      N      n     01c0 111 0 == '7' * 4 == 16 -> 0x...........4.... -> 4
//                        ^^^
//
//  Encoding is technically only defined for the ten letters above.  While
//  you can give it any input, it will happily encode [pq] to 0, [brs] to 1,
//  [deu] to 2, [fvw] to 3 and the rest of the letters to 4 (unless I
//  squinted at the ASCII table wrong).
//
//  Likewise, decoding is technically only defined for 0, 1, 2, 3 and 4.  It'll
//  not-crash for any input, but they're not correct.
//
//  decode2bitBase() is the simpler of the two functions, and is similar to
//  encode2bitBase().  It needs to shift the large constant right by some
//  multiple of 8 so it can then extract the desired 8-bit result.  The shift
//  amount is just the input 'base' (0, 1, 2, 3 or 4).  We mask it to the
//  lower three bits to keep it a valid shift amount after multiplying by 8.
//
//  encode2bitBase() extracts the highlighted three bits from the input,
//  shifts them left one spot -- equivalent to multiplying the decimal
//  equivalent of those three bits by four, by combining the right shift by
//  one and then the left shift by two into one shift -- and uses that to
//  shift the big constant right to move each (four-bit wide) 'encoded'
//  result into the lower four bits.
//
//  The 'simpler' commented out version of encode2bitBase uses a 'compound
//  literal' which is NOT valid C++.  Both GCC and Clang seem to support it,
//  but complain when -pedantic is enabled.  decode2bitBase(), however, is
//  valid C++.
//
//line uint8 encode2bitBase(char  base) { return((uint8 [8]){0, 1, 2, 3, 4, 4, 4, 4}[base >> 1 & 0x07]); }
//line char  decode2bitBase(uint8 base) { return("ACTGNNNN"[base & 0x07]); }
inline uint8 encode2bitBase(char  base) { return (0x0000000044443210llu >> ((base & 0x0e) << 1)) & 0x0f; }
inline char  decode2bitBase(uint8 base) { return (0x4e4e4e4e47544341llu >> ((base & 0x07) << 3)) & 0xff; }

uint32 encode2bitSequence(uint8 *&chunk, char const *seq, uint32 seqLen);
uint32 encode3bitSequence(uint8 *&chunk, char const *seq, uint32 seqLen);
uint32 encode8bitSequence(uint8 *&chunk, char const *seq, uint32 seqLen);

void   decode2bitSequence(uint8 const *chunk, uint32 chunkLen, char *seq, uint32 seqLen);
void   decode3bitSequence(uint8 const *chunk, uint32 chunkLen, char *seq, uint32 seqLen);
void   decode8bitSequence(uint8 const *chunk, uint32 chunkLen, char *seq, uint32 seqLen);

}  //  merylutil::sequence::v1

#endif  //  MERYLUTIL_SEQUENCE_V1_H
