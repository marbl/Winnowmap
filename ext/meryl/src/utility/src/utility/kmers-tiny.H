
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

#ifndef MERYL_UTIL_KMER_KMER_TINY_H
#define MERYL_UTIL_KMER_KMER_TINY_H

#ifndef MERYL_UTIL_KMER_H
#error "include kmers.H, not this."
#endif

//  Definition of a 'small' kmer.

typedef uint128    kmdata;   //  128 bits of kmer data
typedef uint32     kmpref;   //   32 bits of kmer prefix == 6 bits file prefix, 6 (default) suffix prefix
typedef uint32     kmvalu;   //   32 bits of kmer value
typedef uint64     kmcolo;   //   64 bits of kmer color

constexpr kmvalu   kmvalumax = uint32max;
constexpr kmcolo   kmcolomax = uint64max;


class  kmerTiny {
public:
  kmerTiny() {
    _mer = 0;
  };

  ~kmerTiny() {
  };

  static
  void        setSize(uint32 ms, bool beVerbose=false) {
    _merSize    = ms;

    _fullMask   = 0;
    _fullMask   = ~_fullMask;
    _fullMask >>= 8 * sizeof(kmdata) - (ms * 2);

    _leftMask   = 0;
    _leftMask   = ~_leftMask;
    _leftMask >>= 8 * sizeof(kmdata) - (ms * 2 - 2);

    _leftShift  = ((2 * ms - 2) % (8 * sizeof(kmdata)));

    if (beVerbose)
      fprintf(stderr, "Set global kmer size to " F_U32 " (fullMask=0x%s leftMask=0x%s leftShift=" F_U32 ")\n",
              _merSize, toHex(_fullMask), toHex(_leftMask), _leftShift);
  };

  static
  uint32      merSize(void) { return(_merSize); };

  //  Push an ASCII base onto the mer, shifting the mer to the right or left
  //  to make space for the new base.  Unlike the 'standard' two-bit encoding,
  //  these encode bases as A=00, C=01, G=11, T=10.
  //
  //       +---------+-- upper/lower case bit
  //       |         |
  //    A 1000001 a 1100001 == 00
  //    C 1000011 c 1100011 == 01
  //    G 1000111 g 1100111 == 11
  //    T 1010100 t 1110100 == 10
  //                    ||
  //                    ++-- bits used for 2-bit encoding
  //
  void        addR(kmdata base)       { _mer  = (((_mer << 2) & _fullMask) | (((base >> 1) & 0x03llu)          )              );  };
  void        addL(kmdata base)       { _mer  = (((_mer >> 2) & _leftMask) | (((base >> 1) & 0x03llu) ^ 0x02llu) << _leftShift);  };

  //  Reverse-complementation of a kmer involves complementing the bases in
  //  the mer, revesing the order of all the bases, then aligning the bases
  //  to the low-order bits of the word.
  //
  kmdata      reverseComplement(kmdata mer) const {

    //  Complement the bases

    mer ^= build_uint128(0xaaaaaaaaaaaaaaaallu, 0xaaaaaaaaaaaaaaaallu);

    //  Reverse the mer

    mer = ((mer >>  2) & build_uint128(0x3333333333333333llu, 0x3333333333333333llu)) | ((mer <<  2) & build_uint128(0xccccccccccccccccllu, 0xccccccccccccccccllu));
    mer = ((mer >>  4) & build_uint128(0x0f0f0f0f0f0f0f0fllu, 0x0f0f0f0f0f0f0f0fllu)) | ((mer <<  4) & build_uint128(0xf0f0f0f0f0f0f0f0llu, 0xf0f0f0f0f0f0f0f0llu));
    mer = ((mer >>  8) & build_uint128(0x00ff00ff00ff00ffllu, 0x00ff00ff00ff00ffllu)) | ((mer <<  8) & build_uint128(0xff00ff00ff00ff00llu, 0xff00ff00ff00ff00llu));
    mer = ((mer >> 16) & build_uint128(0x0000ffff0000ffffllu, 0x0000ffff0000ffffllu)) | ((mer << 16) & build_uint128(0xffff0000ffff0000llu, 0xffff0000ffff0000llu));
    mer = ((mer >> 32) & build_uint128(0x00000000ffffffffllu, 0x00000000ffffffffllu)) | ((mer << 32) & build_uint128(0xffffffff00000000llu, 0xffffffff00000000llu));
    mer = ((mer >> 64) & build_uint128(0x0000000000000000llu, 0xffffffffffffffffllu)) | ((mer << 64) & build_uint128(0xffffffffffffffffllu, 0x0000000000000000llu));

    //  Shift and mask out the bases not in the mer

    mer >>= 8 * sizeof(kmdata) - _merSize * 2;
    mer  &= _fullMask;

    return(mer);
  };

  kmerTiny   &reverseComplement(void) {
    _mer = reverseComplement(_mer);
    return(*this);
  };

public:
  bool        operator!=(kmerTiny const &r) const { return(_mer != r._mer); };
  bool        operator==(kmerTiny const &r) const { return(_mer == r._mer); };
  bool        operator< (kmerTiny const &r) const { return(_mer <  r._mer); };
  bool        operator> (kmerTiny const &r) const { return(_mer >  r._mer); };
  bool        operator<=(kmerTiny const &r) const { return(_mer <= r._mer); };
  bool        operator>=(kmerTiny const &r) const { return(_mer >= r._mer); };

  bool        isFirst(void)                 const { return(_mer == 0);         };
  bool        isLast(void)                  const { return(_mer == _fullMask); };

  bool        isCanonical(void)             const { return(_mer <= reverseComplement(_mer));  };
  bool        isPalindrome(void)            const { return(_mer == reverseComplement(_mer));  };

  kmerTiny   &operator++()                        {                           _mer++;  return(*this);  };
  kmerTiny    operator++(int)                     { kmerTiny before = *this;  _mer++;  return(before); };

  kmerTiny   &operator--()                        {                           _mer--;  return(*this);  };
  kmerTiny    operator--(int)                     { kmerTiny before = *this;  _mer--;  return(before); };

public:
  char    *toString(char *str) const {
    for (uint32 ii=0; ii<_merSize; ii++) {
      uint32  bb = (((_mer >> (2 * ii)) & 0x03) << 1);
      str[_merSize-ii-1] = (bb == 0x04) ? ('T') : ('A' + bb);
    }
    str[_merSize] = 0;
    return(str);
  };

  void     recanonicalizeACGTorder(void) {
    kmdata  fmer = _mer;
    kmdata  rmer = reverseComplement(_mer);
    kmdata  mask = _mer;

    mask >>= 1;
    mask  &= build_uint128(0x5555555555555555llu, 0x5555555555555555llu);

    fmer ^= mask;      //  Convert from ACTG ordering to ACGT ordering.
    rmer ^= mask;

    if (fmer < rmer)   //  Recompute canonical based on ACGT ordering.
      _mer = fmer;
    else
      _mer = rmer;

    _mer ^= mask;      //  Convert back to ACTG ordering for printing.
  };

  operator kmdata () const {
    return(_mer);
  };

  operator uint64 () const = delete;   //  Explicitly fail of someone tries to convert us to an integer
  operator  int64 () const = delete;   //  instead of to a kmdata.  Without these, a cast to, say, uint64
  operator uint32 () const = delete;   //  would be first convert to kmdata (uint128) then down to uint64.
  operator  int32 () const = delete;   //  With these, you'll either get a compile-time error (because
  operator uint16 () const = delete;   //  these are private) or link time error (because they're not
  operator  int16 () const = delete;   //  defined.

  void     setPrefixSuffix(kmpref prefix, kmdata suffix, uint32 width) {
    _mer   = prefix;
    _mer <<= width;
    _mer  |= suffix;
  };

private:
  void     operator>>=(uint32 x)  { _mer >>= x; };
  void     operator<<=(uint32 x)  { _mer <<= x; };


private:
public:
  kmdata         _mer;

  static uint32  _merSize;     //  number of bases in this mer

  static kmdata  _fullMask;    //  mask to ensure kmer has exactly _merSize bases in it

  static kmdata  _leftMask;    //  mask out the left-most base.
  static uint32  _leftShift;   //  how far to shift a base to append to the left of the kmer
};


typedef kmerTiny kmer;


#endif  //  MERYL_UTIL_KMER_KMER_TINY_H
