
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

#ifndef MERYL_UTIL_KMER_H
#define MERYL_UTIL_KMER_H

#include "types.H"
#include "system.H"
#include "files.H"
#include "bits.H"

//  merSize 1 NOT supported.  Fails _leftShift.

#undef  SHOW_LOAD

#include "kmers-tiny.H"
#include "kmers-histogram.H"

#include "kmers-iterator.H"

#include "kmers-writer.H"
#include "kmers-reader.H"

#include "kmers-iterator.H"
#include "kmers-lookup.H"


#endif  //  MERYL_UTIL_KMER
