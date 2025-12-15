
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

#ifndef MERYLUTIL_BITS_HEXDUMP_V1_H
#define MERYLUTIL_BITS_HEXDUMP_V1_H

#include "types.H"

//
//  Dump DATAlen bytes from DATA in a hex format.
//  It will print W bytes per line, separated into words of 8 bytes.
//  The end of the line will have the ASCII representation of the data.
//
//  00000000  00 01 02 03 04 05 06 07  08 09 0a 0b 0c 0d 0e 0f  '................'
//

namespace merylutil::inline bits::inline v1 {

void
hexDump(FILE *F,
        void *DATA, uint32 DATAlen, uint32 W=32);

}  //  namespace merylutil::bits::v1

#endif  //  MERYLUTIL_BITS_HEXDUMP_V1_H
