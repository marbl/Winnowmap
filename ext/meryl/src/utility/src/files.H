
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

#ifndef MERYLUTIL_FILES_H
#define MERYLUTIL_FILES_H

//  BSDs don't have this, Linux needs it (or, at least, needed it).
#ifndef O_LARGEFILE
#define O_LARGEFILE 0
#endif

#include "files/accessing-v1.H"
#include "files/reading-v1.H"
#include "files/writing-v1.H"

#include "files/buffered-v1-reading.H"
#include "files/buffered-v1-writing.H"

#include "files/compressed-v1.H"
#include "files/compressed-v1-reading.H"
#include "files/compressed-v1-writing.H"

#include "files/memoryMapped-v1.H"

#include "files/fasta-fastq-v1.H"



//  Search CANU_INSTALL_PATH, MERYL_INSTALL_PATH and PATH for a data file.
char const *findSharedFile(char const *binname, char const *relpath, char const *filename);



#endif  //  MERYLUTIL_FILES_H
