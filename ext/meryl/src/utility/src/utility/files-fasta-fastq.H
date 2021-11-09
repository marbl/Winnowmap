
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

#ifndef FILES_FASTA_FASTQ_H
#define FILES_FASTA_FASTQ_H

//  Write a sequence 's' of length 'sl' to file 'f'.
//
//  For FASTA, the sequence will have line breaks every 'bl' letters.  If
//  'bl' is zero, no line breaks are output.
//
//  FASTQ quality values are either:
//     pre-encoded letters ('!' to '~' inclusive) or
//     unencoded integer values (0 to 93 inclusive).
//
//  Both expect 'h' to be an ident line, WITHOUT the leading '>' or '@' and
//  WITHOUT a trailing newline.  The line can use printf() format
//  declarations, and have the optional arguments passed to it.  Try to keep
//  the length under 4 KB.
//
//  FASTQ qv ident lines are just a single '+'.
//
//  The NUL terminating byte is NOT checked for and NOT used.

void outputFASTA(FILE *f, char const *s,                 uint64 sl, uint64 bl, char const *h, ...);
void outputFASTQ(FILE *f, char const *s, char  const *q, uint64 sl,            char const *h, ...);
void outputFASTQ(FILE *f, char const *s, uint8 const *q, uint64 sl,            char const *h, ...);
void outputFASTQ(FILE *f, char const *s, uint8       qv, uint64 sl,            char const *h, ...);


//  A rather complicated output function to output a sequence
//  as a specific format, or as the format the sequence is.
//
//  If only one of asFASTA and asFASTQ is set, the sequence is output
//  according to those flags:
//    asFASTQ with constant QV if hasQuals is false
//    asFASTQ                  if hasQuals is true
//    asFASTA                  regarless of hasQuals
//
//  If asFASTA is the same asFASTQ, the sequence is output:
//    as FASTA if hasQuals is false
//    as FASTQ if hasQuals is true.
//
//  Cheat sheet:
//      hasQuals  asFASTA  asFASTQ
//      --------  -------  -------
//         false    false    false  ->  fasta
//         false    false     true  ->  fastq with constant qv
//         false     true    false  ->  fasta
//         false     true     true  ->  fasta
//
//          true    false    false  ->  fastq
//          true    false     true  ->  fastq
//          true     true    false  ->  fasta
//          true     true     true  ->  fastq
//
void
outputSequence(FILE        *f,
               char  const *outputName,
               char  const *outputBases,
               uint8 const *outputQuals,  uint32  outputLen,
               bool         hasQuals,
               bool         asFASTA,
               bool         asFASTQ,
               uint8        QV);


#endif  //  FILES_FASTA_FASTQ_H
