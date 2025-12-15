
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

#ifndef MERYLUTIL_MATH_MD5_V1
#define MERYLUTIL_MATH_MD5_V1

#include "types.H"

namespace merylutil::inline math::inline v1 {


//  md5sum  md5;
//  while (data)
//    md5.addBlock(data, dataLen);
//  md5.finalize()
//
class md5sum {
public:
  md5sum();
  md5sum(uint8 const *b, uint32 bLen);
  ~md5sum();

  char const *toString(void)        { return dascii; }
  void        toString(char *str);

  void        toWords(uint64 &hi, uint64 &lo) {
    hi = upper;
    lo = lower;
  }
  uint64      hiWord(void)          { return upper; }
  uint64      loWord(void)          { return lower; }

  void        addBlock(uint8 *b, uint64 bLen);

  void        finalize(void);

private:
  uint8       digest[33] = {0};
  char        dascii[33] = {0};

  uint64      upper = uint64zero;  //0xd41d8cd98f00b204;
  uint64      lower = uint64zero;  //0xe9800998ecf8427e;

  void       *context;
};

}  //  merylutil::files::v1

#endif  //  MERYLUTIL_MATH_MD5_V1
