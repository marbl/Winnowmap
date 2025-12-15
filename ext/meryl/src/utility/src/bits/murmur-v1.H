
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

#ifndef MERYLUTIL_BITS_MURMUR_V1_H
#define MERYLUTIL_BITS_MURMUR_V1_H

#include "bits-v1.H"

////////////////////////////////////////
//
//  Geerate something like a murmur3 hash.  Not the official blessed version;
//  It's the MerylUtilMurMur hash.
//  
namespace merylutil::inline bits::inline v1 {

class mumurmur32 {
public:
  void    init(uint32 s) {           //  Reset to compute a new hash.
    _l = 0;
    _h = s;
  }

  void    add(uint32 w) {            //  Add a single word to the hash.
    _l += 1;                         //
                                     //  Technically, _l should be in bytes.
    _k  = w;                         //
    _k *= ac1;                       //  Mix bits in the word.
    // &= uint32max;                 //
    _k  = rotateBitsLeft(_k, ar1);   //
    _k *= ac2;                       //
    // &= uint32max;                 //
                                     //
    _h ^= _k;                        //  Add the mixed bits to the
    _h  = rotateBitsLeft(_h, ar2);   //  preliminary hash.
    _h += _h << 2;                   //    (aka h = 5*h)
    _h += ac3;                       //
    // &= uint32max;                 //
  }

  uint32  mix(void) {                //  Mix bits to generate the final hash value.
    _h ^= _l;                        //
    _h ^= _h >> ms1;                 //
    _h *= mm1;                       //
    // &= uint32max;                 //
    _h ^= _h >> ms2;                 //  NB: the &'s with uint32max do nothing
    _h *= mm2;                       //  in this implementation; they're there
    // &= uint32max;                 //  for consistency with an implementation
    _h ^= _h >> ms3;                 //  that uses larger _k and _h word sizes.

    return _h;
  }

private:
  uint32   _k = 0;
  uint32   _h = 0;
  uint32   _l = 0;

private:
  static constexpr uint64 ac1 = 0xcc9e2d51;
  static constexpr uint32 ar1 = 15;
  static constexpr uint64 ac2 = 0x1b873593;
  static constexpr uint32 ar2 = 13;
  static constexpr uint64 ac3 = 0xe6546b64;

  static constexpr uint32 ms1 = 16;
  static constexpr uint64 mm1 = 0x85ebca6b;
  static constexpr uint32 ms2 = 13;
  static constexpr uint64 mm2 = 0xc2b2ae35;
  static constexpr uint32 ms3 = 16;
};

inline uint32 mumurmur(uint32 *a, uint32 aLen) {
  mumurmur32 mh;

  mh.init(0xb0f57ee3lu);     //  From smasher
  for (uint32 ii=0; ii<aLen; ii++)
    mh.add(a[ii]);
  return mh.mix();
}

}  //  namespace merylutil::bits::v1



#if 0
//  murmur 3 is defined oly for 32-bit and 128-bit hashes.
//  these constants could possibly be used for a 64-bit version;
//  the 128-bit hash computes two 64-bit hashes and mixes those at the end.

  constexpr uint64 c1 = 0x87c37b91114253d5;
  constexpr uint32 r1 = 31;  //33
  constexpr uint64 c2 = 0x4cf5ad432745937f;
  constexpr uint32 r2 = 27;  //31
  constexpr uint64 c3 = 0x52dce729;  (and 0x38495ab5 for the second 64-bit word)

  constexpr uint32 ms1 = 33;
  constexpr uint64 mm1 = 0xff51afd7ed558ccdllu;
  constexpr uint32 ms2 = 33;
  constexpr uint64 mm2 = 0xc4ceb9fe1a85ec53llu;
  constexpr uint32 ms3 = 33;
#endif


#endif  //  MERYLUTIL_BITS_MURMUR_V1_H
