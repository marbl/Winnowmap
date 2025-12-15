
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

#ifndef MERYLUTIL_KMERS_H
#define MERYLUTIL_KMERS_H

//  merSize 1 NOT supported.  Fails _leftShift.

#undef  SHOW_LOAD

////////////////////////////////////////
//
//  Version 1 is used in meryl (and meryl-lookup), merfin and canu.
//
#include "kmers-v1/kmers-tiny.H"
#include "kmers-v1/kmers-histogram.H"

#include "kmers-v1/kmers-iterator.H"

#include "kmers-v1/kmers-writer.H"
#include "kmers-v1/kmers-reader.H"

#include "kmers-v1/kmers-iterator.H"
#include "kmers-v1/kmers-lookup.H"

////////////////////////////////////////
//
//  Version 2 is used in meryl2.
//
#include "kmers-v2/kmers-tiny.H"
#include "kmers-v2/kmers-histogram.H"

#include "kmers-v2/kmers-iterator.H"

#include "kmers-v2/kmers-writer.H"
#include "kmers-v2/kmers-reader.H"

#include "kmers-v2/kmers-iterator.H"
#include "kmers-v2/kmers-lookup.H"

#endif  //  MERYLUTIL_KMERS
