
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

#include "types.H"
#include "strings.H"
#include "system.H"


bool
printFeature(bool enabled, char const *str, bool first=true) {

  if (enabled == false)
    return(first);

  if (first == false)
    printf(",%s", str);
  else
    printf("%s", str);

  return(false);
}


int
main(int argc, char **argv) {
  merylutil::cpuIdent id;

  printf("CPU: '%s' '%s' Model 0x%02x Family 0x%02x Stepping %u\n",
         id._processorOrigin,
         id._processorName,
         id._modelID,
         id._familyID,
         id._stepping);

  bool f;

  printf("Vector Features: [");
  f = printFeature(id.supportsMMX(),   "MMX");
  f = printFeature(id.supportsSSE(),   "SSE", f);
  f = printFeature(id.supportsSSE2(),  "SSE2", f);
  f = printFeature(id.supportsSSE3(),  "SSE3", f);
  f = printFeature(id.supportsSSSE3(), "SSSE3", f);
  f = printFeature(id.supportsSSE4A(), "SSE4A", f);
  f = printFeature(id.supportsSSE4_1(), "SSE4.1", f);
  f = printFeature(id.supportsSSE4_2(), "SSE4.2", f);
  f = printFeature(id.supportsAVX(), "AVX", f);
  f = printFeature(id.supportsAVX2(), "AVX2", f);
  f = printFeature(id.supportsAMX(), "AMX", f);
  f = printFeature(id.supportsPOPCNT(), "POPCNT", f);
  printf("]\n");

  printf("AVX512 Features: [");
  f = printFeature(id.supportsAVX512(), "AVX512", f);
  f = printFeature(id.supportsAVX512F(), "F");
  f = printFeature(id.supportsAVX512DQ(), "DQ", f);
  f = printFeature(id.supportsAVX512IFMA(), "IFMA", f);
  f = printFeature(id.supportsAVX512PF(), "PF", f);
  f = printFeature(id.supportsAVX512ER(), "ER", f);
  f = printFeature(id.supportsAVX512CD(), "CD", f);
  f = printFeature(id.supportsAVX512BW(), "BW", f);
  f = printFeature(id.supportsAVX512VL(), "VL", f);
  f = printFeature(id.supportsAVX512VBMI(), "VBMI", f);
  f = printFeature(id.supportsAVX512VBMI2(), "VBMI2", f);
  f = printFeature(id.supportsAVX512VNNI(), "VNNI", f);
  f = printFeature(id.supportsAVX512BITALG(), "BITALG", f);
  f = printFeature(id.supportsAVX512VPOPCNTDQ(), "VPOPCNTDQ", f);
  f = printFeature(id.supportsAVX5124VNNIW(), "4VNNIW", f);
  f = printFeature(id.supportsAVX5124FMAPS(), "4FMAPS", f);
  f = printFeature(id.supportsAVX512BF16(), "BF16", f);
  printf("]\n");

  return 0;
}

