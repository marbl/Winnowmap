
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

#ifdef __x86_64__
#include <cpuid.h>
#include <x86intrin.h>
#include <immintrin.h>
#endif


//https://android.googlesource.com/platform/ndk/+/ndk-release-r21/sources/android/cpufeatures/cpu-features.c
#ifdef __aarch64__
//#include <uapi/asm/hwcap.h>
#endif

//https://developer.arm.com/documentation/ddi0406/c/System-Level-Architecture/The-CPUID-Identification-Scheme/The-CPUID-registers/Organization-of-the-CPUID-registers?lang=en

#ifdef __arm__
//#include <asm/hwcap.h>
//#include <uapi/asm/hwcap.h>
#endif

#ifdef __mips__
#endif

#include "cpuIdent-v1.H"

//
//  Load processor flags.  They'll be interpreted in the various has*() functions.
//
void
merylutil::system::v1::cpuIdent::loadProcessorFlags(void) {

#ifdef __x86_64__
  __builtin_cpu_init();

  _maxLeaf    = __get_cpuid_max(0x00000000, nullptr);
  _maxExtLeaf = __get_cpuid_max(0x80000000, nullptr);

  if (isLeafValid(0x01)) __cpuid_count(0x01, 0, _leaf1sub0eax, _leaf1sub0ebx, _leaf1sub0ecx, _leaf1sub0edx);

  if (isLeafValid(0x07)) __cpuid_count(0x07, 0, _leaf7sub0eax, _leaf7sub0ebx, _leaf7sub0ecx, _leaf7sub0edx);
  if (isLeafValid(0x07)) __cpuid_count(0x07, 1, _leaf7sub1eax, _leaf7sub1ebx, _leaf7sub1ecx, _leaf7sub1edx);

  if (isExtendedLeafValid(0x01)) __cpuid_count(0x80000001, 0, _extf1sub0eax, _extf1sub0ebx, _extf1sub0ecx, _extf1sub0edx);

  if (hasXGETBV())
    _xcr0 = _xgetbv(0);   //  Needs clang/gcc option -mxsave
#endif
}


//
//  Get the processor origin and decide if we are Intel or AMD.
//
void
merylutil::system::v1::cpuIdent::decodeProcessorOrigin(void) {
  uint32  *r = (uint32 *)_processorOrigin;

  r[000] = r[001] = r[002] = r[003] = 0;

#ifdef __x86_64__
  if (isLeafValid(0) == false)
    return;

  __cpuid_count(0x00000000, 0, r[003], r[000], r[002], r[001]);  r[003] = 0;

  _isIntel = ((r[000] == 0x756e6547) && (r[001] == 0x49656e69) && (r[002] == 0x6c65746e));
  _isAMD   = ((r[000] == 0x68747541) && (r[001] == 0x69746e65) && (r[002] == 0x444d4163));

  merylutil::trimString(_processorOrigin);
  //printf("'%s' 0x%08x 0x%08x 0x%08x 0x%08x\n", _processorOrigin, r[0], r[1], r[2], r[3]);
#endif
}


//
//  Get the processor name.
//
//  Examples:
//    'AMD Ryzen 7 3700X 8-Core Processor             '  - 
//    'Genuine Intel(R) CPU 0000 @ 2.20GHz'              - Intel E5-2650 v4 ES
//    '        Intel(R) Atom(TM) CPU D2500   @ 1.86GHz'  - Intel D2500 SBC
//    'Intel(R) Core(TM) i7-4870HQ CPU @ 2.50GHz'        - 2016-ish MacBook Pro
//    'Intel(R) Core(TM) i7-9750H CPU @ 2.60GHz'         - 2020-ish MacBook Pro
//
void
merylutil::system::v1::cpuIdent::decodeProcessorName(void) {
  uint32  *r = (uint32 *)_processorName;

  r[000] = r[001] = r[002] = r[003] = 0;
  r[004] = r[005] = r[006] = r[007] = 0;
  r[010] = r[011] = r[012] = r[013] = 0;

#ifdef __x86_64__
  if (isExtendedLeafValid(4) == false)
    return;

  __cpuid_count(0x80000002, 0, r[000], r[001], r[002], r[003]);
  __cpuid_count(0x80000003, 0, r[004], r[005], r[006], r[007]);
  __cpuid_count(0x80000004, 0, r[010], r[011], r[012], r[013]);

  merylutil::trimString(_processorName);
#endif
}


//
//  Decode the proccessor model, vendor and stepping.
//
void
merylutil::system::v1::cpuIdent::decodeProcessorModelVendorStepID(void) {

#ifdef __x86_64__
  if (isLeafValid(1) == false)   return;

  int extFamID   = (_leaf1sub0eax >> 20) & 0xff;
  int extModID   = (_leaf1sub0eax >> 16) & 0x0f;
  int procType   = (_leaf1sub0eax >> 12) & 0x03;
  int familyID   = (_leaf1sub0eax >>  8) & 0x0f;
  int modelID    = (_leaf1sub0eax >>  4) & 0x0f;
  int stepID     = (_leaf1sub0eax >>  0) & 0x0f;

  int logProcID  = (_leaf1sub0ebx >> 24) & 0xff;
  int nLogProc   = (_leaf1sub0ebx >> 16) & 0xff;
  int clFlushLS  = (_leaf1sub0ebx >>  8) & 0xff;
  int brandIndex = (_leaf1sub0ebx >>  0) & 0xff;

  //fprintf(stdout, "logProc %d nLogProc %d clFlushLS %d brandIndex %d\n", logProc, nLogProc, clFlushLS, brandIndex);

  if ((familyID == 6) || (familyID == 15))
    modelID = (extModID << 4) | (modelID);

  if (familyID == 15)
    familyID = extFamID + familyID;

  //printf("CPU: ID=0x%08x Family=0x%02x Model=0x%02x Stepping=0x%02x\n", eax, familyID, modelID, stepID);

  _modelID  = modelID;
  _familyID = familyID;
  _stepping = stepID;
#endif
}


//
//  Decode Intel processor topology.
//
bool
merylutil::system::v1::cpuIdent::decodeIntelCacheTopology(void) {
  uint32  eax, ebx, ecx, edx;

  if (false == _isIntel)           return false;
  if (false == isLeafValid(0x04))  return false;

#ifdef __x86_64__
  for (uint32 c=0; c<32; c++) {
    __cpuid_count(0x04, c, eax, ebx, ecx, edx);

    uint32 cType          =  (eax >>  0) & 0x001f;
    uint32 cLevel         =  (eax >>  5) & 0x0007;
    uint32 cInit          =  (eax >>  8) & 0x0001;
    uint32 cFAssoc        =  (eax >>  9) & 0x0001;
    //nt32 reserved       =  (eax >> 10) & 0x000f;
    uint32 nLogProc       = ((eax >> 14) & 0x0fff) + 1;
    uint32 nPhysProc      = ((eax >> 26) & 0x003f) + 1;

    uint32 cLineSize      = ((ebx >>  0) & 0x0fff) + 1;
    uint32 cPartitions    = ((ebx >> 12) & 0x03ff) + 1;
    uint32 cAssociativity = ((ebx >> 22) & 0x03ff) + 1;

    uint32 nSets          =   ecx + 1;

    uint32 noInvalidate   =  (edx >>  0) & 0x0001;
    uint32 inclusive      =  (edx >>  1) & 0x0002;
    uint32 complexIndex   =  (edx >>  2) & 0x0004;

    if (cType == 0x00)
      break;

    uint32 cSize = cAssociativity * cPartitions * cLineSize * nSets;

    if (cType == 0x01)  fprintf(stderr, "Intel L%u Data     %5uKB %2d-way per %2d logical processors\n", cLevel, cSize >> 10, cAssociativity, nLogProc);
    if (cType == 0x02)  fprintf(stderr, "Intel L%u Instr    %5uKB %2d-way per %2d logical processors\n", cLevel, cSize >> 10, cAssociativity, nLogProc);
    if (cType == 0x03)  fprintf(stderr, "Intel L%u Unified  %5uKB %2d-way per %2d logical processors\n", cLevel, cSize >> 10, cAssociativity, nLogProc);
  }
#endif

  return true;
}


//
//  Decode AMD's newer processor topology.  This is almost the same as Intel, just a few
//  undefined pieces.
//
//  Described on PDF page 671 in AMD64 Architecture Programmer's Manual
//  Volume 3, General Purpose and System Instructions.
//
bool
merylutil::system::v1::cpuIdent::decodeAMDCacheTopology(void) {
  uint32 eax, ebx, ecx, edx;

  if (false == _isAMD)                     return false;
  if (false == isExtendedLeafValid(0x1d))  return false;
  if (false == hasTOPOEXT())               return false;

#ifdef __x86_64__
  for (uint32 c=0; c<32; c++) {
    __cpuid_count(0x8000001d, c, eax, ebx, ecx, edx);

    uint32 cType           =  (eax >>  0) & 0x001f;
    uint32 cLevel          =  (eax >>  5) & 0x0007;
    uint32 cInit           =  (eax >>  8) & 0x0001;
    uint32 cFAssoc         =  (eax >>  9) & 0x0001;
    //nt32 reserved        =  (eax >> 10) & 0x000f;
    uint32 nLogProc        = ((eax >> 14) & 0x0fff) + 1;
    //nt32 nPhysProc       =

    uint32 cLineSize       = ((ebx >>  0) & 0x0fff) + 1;
    uint32 cPartitions     = ((ebx >> 12) & 0x03ff) + 1;
    uint32 cAssociativity  = ((ebx >> 22) & 0x03ff) + 1;

    uint32 nSets           =   ecx + 1;

    uint32 noInvalidate    =  (edx >>  0) & 0x0001;
    uint32 inclusive       =  (edx >>  1) & 0x0002;
    //nt32 complexIndex    =

    if (cType == 0x00)
      break;

    uint32 cSize = cAssociativity * cPartitions * cLineSize * nSets;

    if (cType == 0x01)  fprintf(stderr, "AMD   L%u Data     %5uKB %2d-way per %2d logical processors\n", cLevel, cSize >> 10, cAssociativity, nLogProc);
    if (cType == 0x02)  fprintf(stderr, "AMD   L%u Instr    %5uKB %2d-way per %2d logical processors\n", cLevel, cSize >> 10, cAssociativity, nLogProc);
    if (cType == 0x03)  fprintf(stderr, "AMD   L%u Unified  %5uKB %2d-way per %2d logical processors\n", cLevel, cSize >> 10, cAssociativity, nLogProc);
  }
#endif

  return true;
}


bool
merylutil::system::v1::cpuIdent::decodeAMDLegacyCacheTopology(void) {
  uint32 eax, ebx, ecx, edx;

  if (false == _isAMD)             return false;
  if (false == isLeafValid(0x05))  return false;
  if (false == isLeafValid(0x06))  return false;

#ifdef __x86_64__
  uint32  assocTable[16] = {0, 1, 2, 0, 4, 0, 8, 0, 16, 0, 32, 48, 64, 96, 128, uint32max};

  __cpuid_count(0x80000005, 0, eax, ebx, ecx, edx);

  uint32  cL1DataLineSize      =  (ecx >>  0) & 0x00ff;
  uint32  cL1DataLinePerTag    =  (ecx >>  8) & 0x00ff;
  uint32  cL1DataAssociativity =  (ecx >> 16) & 0x00ff;
  uint32  cL1DataSizeKB        =  (ecx >> 24) & 0x00ff;

  uint32  cL1InstLineSize      =  (edx >>  0) & 0x00ff;
  uint32  cL1InstLinePerTag    =  (edx >>  8) & 0x00ff;
  uint32  cL1InstAssociativity =  (edx >> 16) & 0x00ff;
  uint32  cL1InstSizeKB        =  (edx >> 24) & 0x00ff;

  __cpuid_count(0x80000006, 0, eax, ebx, ecx, edx);

  uint32  cL2LineSize      =  (ecx >>  0) & 0x00ff;
  uint32  cL2LinePerTag    =  (ecx >>  8) & 0x000f;
  uint32  cL2Associativity =  (ecx >> 12) & 0x000f;
  uint32  cL2SizeKB        =  (ecx >> 16) & 0xffff;

  cL2Associativity = assocTable[cL2Associativity];

  uint32  cL3LineSize      =  (edx >>  0) & 0x00ff;
  uint32  cL3LinePerTag    =  (edx >>  8) & 0x000f;
  uint32  cL3Associativity =  (edx >> 12) & 0x000f;
  uint32  cL3SizeKB        = ((edx >> 18) & 0x3fff) * 512;

  assert(cL3Associativity != 9);  //  New method required!

  cL3Associativity = assocTable[cL3Associativity];

  fprintf(stderr, "AMDold L1 Data     %5uKB %3u B/line %d-way\n", cL1DataSizeKB, cL1DataLineSize, cL1DataAssociativity);
  fprintf(stderr, "AMDold L1 Instr    %5uKB %3u B/line %d-way\n", cL1InstSizeKB, cL1InstLineSize, cL1InstAssociativity);
  fprintf(stderr, "AMDold L2 Unified  %5uKB %3u B/line %d-way\n", cL2SizeKB, cL2LineSize, cL2Associativity);
  fprintf(stderr, "AMDold L3 Unified  %5uKB %3u B/line %d-way\n", cL3SizeKB, cL3LineSize, cL3Associativity);
#endif

  return true;
}

