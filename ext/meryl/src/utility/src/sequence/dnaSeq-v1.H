
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

#ifndef MERYLUTIL_SEQUENCE_V1_DNASEQ_H
#define MERYLUTIL_SEQUENCE_V1_DNASEQ_H

#include "types.H"

//
//  A sequence loaded from disk.  It should be treated as a read-only object.
//
//  header() returns the entire unmodified sequence header line,
//  in particular, no spaces are removed (see ident() and flags()).
//
//  ident() returns the first word of the sequence header line.
//   - Spaces before the first word are always stripped.
//
//  flags() returns the rest of the line, or an empty line if there are no flags.
//   - If the first whitespace after the ident word is a TAB,
//     then the rest of the header line is returned as the flags:
//       header() == \s* + ident() + '\t' + flags()
//   - Otherwise, whitespace between the ident and the first flag
//     word and all spaces at the end of the line are removed.
//       header() == \s* + ident() + \s* + flags() + \s*
//
//  ident() and flags() could support modification (they're independent stings)
//  but it would be tedious to support length changes and nobody does it anyway.
//
//  bases() and quals() could support modifications, as long as the length of
//  the string doesn't change.  Only canu needed to do that and it was worked
//  around.
//
//  If quality values are not available (e.g., FASTA) then all values are set
//  to zero.
//
//  The copy functions will copy bases (and qualities) from bgn to end, but
//  not including the base at end -- that is, normal C-style semantics.  The
//  output will be NUL-terminated, unless explicitly told not to.  Returns
//  false if bgn or end are out of range or inconsistent.
//

namespace merylutil::inline sequence::inline v1 {

class dnaSeq {
public:
  dnaSeq();
  ~dnaSeq();

  char  const  *header(void)       { return _headr;  };
  char  const  *ident(void)        { return _ident;  };
  char  const  *flags(void)        { return _flags;  };

  char  const  *bases(void)        { return _seq;    };
  uint8 const  *quals(void)        { return _qlt;    };

  uint64        length(void)       { return _seqLen; };

private:
  void          releaseIdent(void);    //  Release ident memory.
  void          releaseBases(void);    //  Release seq/qlt memory.

public:
  bool          copy(char  *bout,
                     uint32 bgn, uint32 end, bool terminate = true);

  bool          copy(char  *bout,
                     uint8 *qout,
                     uint32 bgn, uint32 end, bool terminate = true);

  bool          wasError(void)     { return (_error & 0x01) == 0x01; };
  bool          wasReSync(void)    { return (_error & 0x02) == 0x02; };

private:
  void          findNameAndFlags(void);

private:
  char         *_headr  = nullptr;   uint32  _headrMax = 0;
  char         *_ident  = nullptr;   uint32  _identMax = 0;
  char         *_flags  = nullptr;   uint32  _flagsMax = 0;

  char         *_seq    = nullptr;
  uint8        *_qlt    = nullptr;
  uint64        _seqMax = 0;         //  Space allocated.
  uint64        _seqLen = 0;         //  Actual length.

  uint32        _error  = 0;

  friend class dnaSeqFile;
};

}  //  namespace merylutil::sequence::v1

#endif  //  MERYLUTIL_SEQUENCE_V1_DNASEQ_H
