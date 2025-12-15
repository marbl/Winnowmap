
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

#ifndef MERYL_H
#define MERYL_H

#include <stack>
#include <vector>

#include "system.H"
#include "strings.H"
#include "bits.H"
#include "kmers.H"
#include "sequence.H"


using namespace merylutil;
using namespace merylutil::kmers::v2;

#define MERYLINCLUDE

#include "merylGlobals.H"
#include "merylInput.H"

#include "merylAssign.H"
#include "merylSelector.H"

#include "merylCommandBuilder.H"

#include "merylOp.H"
#include "merylOpCounting.H"
#include "merylOpTemplate.H"
#include "merylOpCompute.H"

#include "merylCountArray.H"

#undef  MERYLINCLUDE

#endif  //  MERYL_H
