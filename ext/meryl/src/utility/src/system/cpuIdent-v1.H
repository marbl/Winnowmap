
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

#ifndef MERYLUTIL_SYSTEM_CPUID_H
#define MERYLUTIL_SYSTEM_CPUID_H

#include "types.H"
#include "strings.H"

namespace merylutil::inline system::inline v1 {

class cpuIdent {
public:
  cpuIdent() {
    loadProcessorFlags();

    decodeProcessorOrigin();
    decodeProcessorName();

    decodeProcessorModelVendorStepID();

    bool validCacheTopology = (decodeIntelCacheTopology() || decodeAMDCacheTopology() || decodeAMDLegacyCacheTopology());
  }
  ~cpuIdent() {
  }

public:

  bool supportsMMX(void)             { return hasMMX();                                  }
  bool supportsSSE(void)             { return hasSSE()    && hasXSAVE() && hasOSXSAVE(); }   //  not sure if XSAVE needed or just fxsr
  bool supportsSSE2(void)            { return hasSSE2()   && hasXSAVE() && hasOSXSAVE(); }
  bool supportsSSE3(void)            { return hasSSE3()   && hasXSAVE() && hasOSXSAVE(); }
  bool supportsSSSE3(void)           { return hasSSSE3()  && hasXSAVE() && hasOSXSAVE(); }
  bool supportsSSE4A(void)           { return hasSSE4A()  && hasXSAVE() && hasOSXSAVE(); }
  bool supportsSSE4_1(void)          { return hasSSE4_1() && hasXSAVE() && hasOSXSAVE(); }
  bool supportsSSE4_2(void)          { return hasSSE4_2() && hasXSAVE() && hasOSXSAVE(); }

  //  AVX needs XSAVE and OSXSAVE, and 6 in xcr0 (== SSE and AVX)
  //  AVX512 needs 'avx512save' == 0xe0 in xcr0 (== OPMASK, ZMMHI256, HI16ZMM

  bool supportsAVX(void)             { return supportsSSE() && hasAVX();  }
  bool supportsAVX2(void)            { return supportsAVX() && hasAVX2(); }

  bool supportsVPCLMULQDQ(void)      { return supportsAVX() && hasVPCLMULQDQ(); }

  bool supportsAVX512(void)          { return supportsSSE() && hasAVX() && hasOPMASK() && hasZMMHI256() && hasHI16ZMM(); }

  bool supportsAVX512F(void)         { return supportsAVX512() && hasAVX512F();         }
  bool supportsAVX512DQ(void)        { return supportsAVX512() && hasAVX512DQ();        }
  bool supportsAVX512IFMA(void)      { return supportsAVX512() && hasAVX512IFMA();      }
  bool supportsAVX512PF(void)        { return supportsAVX512() && hasAVX512PF();        }
  bool supportsAVX512ER(void)        { return supportsAVX512() && hasAVX512ER();        }
  bool supportsAVX512CD(void)        { return supportsAVX512() && hasAVX512CD();        }
  bool supportsAVX512BW(void)        { return supportsAVX512() && hasAVX512BW();        }
  bool supportsAVX512VL(void)        { return supportsAVX512() && hasAVX512VL();        }
  bool supportsAVX512VBMI(void)      { return supportsAVX512() && hasAVX512VBMI();      }
  bool supportsAVX512VBMI2(void)     { return supportsAVX512() && hasAVX512VBMI2();     }
  bool supportsAVX512VNNI(void)      { return supportsAVX512() && hasAVX512VNNI();      }
  bool supportsAVX512BITALG(void)    { return supportsAVX512() && hasAVX512BITALG();    }
  bool supportsAVX512VPOPCNTDQ(void) { return supportsAVX512() && hasAVX512VPOPCNTDQ(); }
  bool supportsAVX5124VNNIW(void)    { return supportsAVX512() && hasAVX5124VNNIW();    }
  bool supportsAVX5124FMAPS(void)    { return supportsAVX512() && hasAVX5124FMAPS();    }
  bool supportsAVX512BF16(void)      { return supportsAVX512() && hasAVX512BF16();      }

  bool supportsAMX(void)             { return hasTILECFG() && hasTILEDATA(); }

  bool supportsPOPCNT(void)          { return hasPOPCNT(); }

private:

  //  Extended Control Register 0

  bool hasX87(void)             { return ((_xcr0 & 0x00000001) != 0); }
  //ol hasSSE(void)             { return ((_xcr0 & 0x00000002) != 0); }
  //ol hasAVX(void)             { return ((_xcr0 & 0x00000004) != 0); }
  bool hasBNDREG(void)          { return ((_xcr0 & 0x00000008) != 0); }
  bool hasBNDCSR(void)          { return ((_xcr0 & 0x00000010) != 0); }
  bool hasOPMASK(void)          { return ((_xcr0 & 0x00000020) != 0); }
  bool hasZMMHI256 (void)       { return ((_xcr0 & 0x00000040) != 0); }
  bool hasHI16ZMM(void)         { return ((_xcr0 & 0x00000080) != 0); }
  //ol hasRESERVED(void)        { return ((_xcr0 & 0x00000100) != 0); }
  bool hasPKRU(void)            { return ((_xcr0 & 0x00000200) != 0); }
  //ol hasRESERVED(void)        { return ((_xcr0 & 0x00000400) != 0); }
  bool hasTILECFG(void)         { return ((_xcr0 & 0x00010000) != 0); }
  bool hasTILEDATA(void)        { return ((_xcr0 & 0x00020000) != 0); }

  //  Leaf 1 ECX

  bool hasSSE3(void)            { return ((_leaf1sub0ecx & 0x00000001) != 0); }
  bool hasPCLMULQDQ(void)       { return ((_leaf1sub0ecx & 0x00000002) != 0); }   //  PCLMUL
  bool hasDTES64(void)          { return ((_leaf1sub0ecx & 0x00000004) != 0); }
  bool hasMON(void)             { return ((_leaf1sub0ecx & 0x00000008) != 0); }   //  MONITOR
  bool hasDSCPL(void)           { return ((_leaf1sub0ecx & 0x00000010) != 0); }
  bool hasVMX(void)             { return ((_leaf1sub0ecx & 0x00000020) != 0); }   //  LZCNT
  bool hasSMX(void)             { return ((_leaf1sub0ecx & 0x00000040) != 0); }
  bool hasEIST(void)            { return ((_leaf1sub0ecx & 0x00000080) != 0); }
  bool hasTM2(void)             { return ((_leaf1sub0ecx & 0x00000100) != 0); }
  bool hasSSSE3(void)           { return ((_leaf1sub0ecx & 0x00000200) != 0); }
  bool hasCNXTID(void)          { return ((_leaf1sub0ecx & 0x00000400) != 0); }
  bool hasSDBG(void)            { return ((_leaf1sub0ecx & 0x00000800) != 0); }   //  Silicon Debug interface
  bool hasFMA(void)             { return ((_leaf1sub0ecx & 0x00001000) != 0); }
  bool hasCX16(void)            { return ((_leaf1sub0ecx & 0x00002000) != 0); }   //  CMPXCHG16B
  bool hasXTPR(void)            { return ((_leaf1sub0ecx & 0x00004000) != 0); }
  bool hasPDCM(void)            { return ((_leaf1sub0ecx & 0x00008000) != 0); }
  //ol hasRESERVED(void)        { return ((_leaf1sub0ecx & 0x00010000) != 0); }   //  reserved
  bool hasPCID(void)            { return ((_leaf1sub0ecx & 0x00020000) != 0); }
  bool hasDCA(void)             { return ((_leaf1sub0ecx & 0x00040000) != 0); }
  bool hasSSE4_1(void)          { return ((_leaf1sub0ecx & 0x00080000) != 0); }
  bool hasSSE4_2(void)          { return ((_leaf1sub0ecx & 0x00100000) != 0); }
  bool hasX2APIC(void)          { return ((_leaf1sub0ecx & 0x00200000) != 0); }
  bool hasMOVBE(void)           { return ((_leaf1sub0ecx & 0x00400000) != 0); }
  bool hasPOPCNT(void)          { return ((_leaf1sub0ecx & 0x00800000) != 0); }
  bool hasTSCDeadline(void)     { return ((_leaf1sub0ecx & 0x01000000) != 0); }
  bool hasAESNI(void)           { return ((_leaf1sub0ecx & 0x02000000) != 0); }   //  AES
  bool hasXSAVE(void)           { return ((_leaf1sub0ecx & 0x04000000) != 0); }
  bool hasOSXSAVE(void)         { return ((_leaf1sub0ecx & 0x08000000) != 0); }
  bool hasAVX(void)             { return ((_leaf1sub0ecx & 0x10000000) != 0); }
  bool hasF16C(void)            { return ((_leaf1sub0ecx & 0x20000000) != 0); }
  bool hasRDRAND(void)          { return ((_leaf1sub0ecx & 0x40000000) != 0); }   //  RDRND
  bool hasHYPERVISOR(void)      { return ((_leaf1sub0ecx & 0x80000000) != 0); }   //  Running under a hypervisor

  bool hasXGETBV(void)          { return hasOSXSAVE() && hasAVX(); }

  //  Leaf 1 EDX

  bool hasFPU()                 { return ((_leaf1sub0edx & 0x00000001) != 0); }
  bool hasVME()                 { return ((_leaf1sub0edx & 0x00000002) != 0); }
  bool hasDE()                  { return ((_leaf1sub0edx & 0x00000004) != 0); }
  bool hasPSE()                 { return ((_leaf1sub0edx & 0x00000008) != 0); }
  bool hasTSC()                 { return ((_leaf1sub0edx & 0x00000010) != 0); }
  bool hasMSR()                 { return ((_leaf1sub0edx & 0x00000020) != 0); }
  bool hasPAE()                 { return ((_leaf1sub0edx & 0x00000040) != 0); }
  bool hasMCE()                 { return ((_leaf1sub0edx & 0x00000080) != 0); }
  bool hasCX8()                 { return ((_leaf1sub0edx & 0x00000100) != 0); }   //  CMPXCHG8B
  bool hasAPIC()                { return ((_leaf1sub0edx & 0x00000200) != 0); }
  //ol hasRESERVED()            { return ((_leaf1sub0edx & 0x00000400) != 0); }   //  reserved
  bool hasSEP()                 { return ((_leaf1sub0edx & 0x00000800) != 0); }
  bool hasMTRR()                { return ((_leaf1sub0edx & 0x00001000) != 0); }
  bool hasPGE()                 { return ((_leaf1sub0edx & 0x00002000) != 0); }
  bool hasMCA()                 { return ((_leaf1sub0edx & 0x00004000) != 0); }
  bool hasCMOV()                { return ((_leaf1sub0edx & 0x00008000) != 0); }
  bool hasPAT()                 { return ((_leaf1sub0edx & 0x00010000) != 0); }
  bool hasPSE36()               { return ((_leaf1sub0edx & 0x00020000) != 0); }
  bool hasPSN()                 { return ((_leaf1sub0edx & 0x00040000) != 0); }
  bool hasCLFLUSH()             { return ((_leaf1sub0edx & 0x00080000) != 0); }   //  CLFSH
  //ol hasRESERVED()            { return ((_leaf1sub0edx & 0x00100000) != 0); }   //  reserved
  bool hasDS()                  { return ((_leaf1sub0edx & 0x00200000) != 0); }
  bool hasACPI()                { return ((_leaf1sub0edx & 0x00400000) != 0); }
  bool hasMMX()                 { return ((_leaf1sub0edx & 0x00800000) != 0); }
  bool hasFXSR()                { return ((_leaf1sub0edx & 0x01000000) != 0); }   //  FXSAVE
  bool hasSSE()                 { return ((_leaf1sub0edx & 0x02000000) != 0); }
  bool hasSSE2()                { return ((_leaf1sub0edx & 0x04000000) != 0); }
  bool hasSS()                  { return ((_leaf1sub0edx & 0x08000000) != 0); }
  bool hasHTT()                 { return ((_leaf1sub0edx & 0x10000000) != 0); }
  bool hasTM()                  { return ((_leaf1sub0edx & 0x20000000) != 0); }
  bool hasIA64()                { return ((_leaf1sub0edx & 0x40000000) != 0); }   //  IA64 processor running x86 emulation
  bool hasPBE()                 { return ((_leaf1sub0edx & 0x40000000) != 0); }

  //  Leaf 7 SubLeaf 0 EBX

  bool hasFSGSBASE(void)        { return ((_leaf7sub0ebx & 0x00000001) != 0); }
  bool hasSGX(void)             { return ((_leaf7sub0ebx & 0x00000004) != 0); }
  bool hasBMI1(void)            { return ((_leaf7sub0ebx & 0x00000008) != 0); }   //  BMI
  bool hasHLE(void)             { return ((_leaf7sub0ebx & 0x00000010) != 0); }
  bool hasAVX2(void)            { return ((_leaf7sub0ebx & 0x00000020) != 0); }
  bool hasSMEP(void)            { return ((_leaf7sub0ebx & 0x00000080) != 0); }
  bool hasBMI2(void)            { return ((_leaf7sub0ebx & 0x00000100) != 0); }
  bool hasENH_MOVSB(void)       { return ((_leaf7sub0ebx & 0x00000200) != 0); }
  bool hasINVPCID(void)         { return ((_leaf7sub0ebx & 0x00000400) != 0); }
  bool hasRTM(void)             { return ((_leaf7sub0ebx & 0x00000800) != 0); }
  bool hasMPX(void)             { return ((_leaf7sub0ebx & 0x00004000) != 0); }
  bool hasAVX512F(void)         { return ((_leaf7sub0ebx & 0x00010000) != 0); }
  bool hasAVX512DQ(void)        { return ((_leaf7sub0ebx & 0x00020000) != 0); }
  bool hasRDSEED(void)          { return ((_leaf7sub0ebx & 0x00040000) != 0); }
  bool hasADX(void)             { return ((_leaf7sub0ebx & 0x00080000) != 0); }
  bool hasAVX512IFMA(void)      { return ((_leaf7sub0ebx & 0x00200000) != 0); }
  bool hasCLFLUSHOPT(void)      { return ((_leaf7sub0ebx & 0x00800000) != 0); }
  bool hasCLWB(void)            { return ((_leaf7sub0ebx & 0x01000000) != 0); }
  bool hasAVX512PF(void)        { return ((_leaf7sub0ebx & 0x04000000) != 0); }
  bool hasAVX512ER(void)        { return ((_leaf7sub0ebx & 0x08000000) != 0); }
  bool hasAVX512CD(void)        { return ((_leaf7sub0ebx & 0x10000000) != 0); }
  bool hasSHA(void)             { return ((_leaf7sub0ebx & 0x20000000) != 0); }
  bool hasAVX512BW(void)        { return ((_leaf7sub0ebx & 0x40000000) != 0); }
  bool hasAVX512VL(void)        { return ((_leaf7sub0ebx & 0x80000000) != 0); }

  //  Leaf 7 SubLeaf 0 ECX

  bool hasPREFTCHWT1(void)      { return ((_leaf7sub0ecx & 0x00000001) != 0); }
  bool hasAVX512VBMI(void)      { return ((_leaf7sub0ecx & 0x00000002) != 0); }
  bool hasUMIP(void)            { return ((_leaf7sub0ecx & 0x00000004) != 0); }   //  PKU
  //ol hasRESERVED()            { return ((_leaf7sub0ecx & 0x00000008) != 0); }
  bool hasOSPKE(void)           { return ((_leaf7sub0ecx & 0x00000010) != 0); }
  bool hasWAITPKG(void)         { return ((_leaf7sub0ecx & 0x00000020) != 0); }
  bool hasAVX512VBMI2(void)     { return ((_leaf7sub0ecx & 0x00000040) != 0); }
  bool hasSHSTK(void)           { return ((_leaf7sub0ecx & 0x00000080) != 0); }
  bool hasGFNI(void)            { return ((_leaf7sub0ecx & 0x00000100) != 0); }
  bool hasVAES(void)            { return ((_leaf7sub0ecx & 0x00000200) != 0); }
  bool hasVPCLMULQDQ(void)      { return ((_leaf7sub0ecx & 0x00000400) != 0); }
  bool hasAVX512VNNI(void)      { return ((_leaf7sub0ecx & 0x00000800) != 0); }
  bool hasAVX512BITALG(void)    { return ((_leaf7sub0ecx & 0x00001000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0ecx & 0x00002000) != 0); }
  bool hasAVX512VPOPCNTDQ(void) { return ((_leaf7sub0ecx & 0x00004000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0ecx & 0x00008000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0ecx & 0x00010000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0ecx & 0x00020000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0ecx & 0x00040000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0ecx & 0x00080000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0ecx & 0x00100000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0ecx & 0x00200000) != 0); }
  bool hasRDPID(void)           { return ((_leaf7sub0ecx & 0x00400000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0ecx & 0x00800000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0ecx & 0x01000000) != 0); }
  bool hasCLDEMOTE(void)        { return ((_leaf7sub0ecx & 0x02000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0ecx & 0x04000000) != 0); }
  bool hasMOVDIRI(void)         { return ((_leaf7sub0ecx & 0x08000000) != 0); }
  bool hasMOVDIR64B(void)       { return ((_leaf7sub0ecx & 0x10000000) != 0); }
  bool hasENQCMD(void)          { return ((_leaf7sub0ecx & 0x20000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0ecx & 0x40000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0ecx & 0x80000000) != 0); }

  //  Leaf 7 SubLeaf 0 EDX

  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x00000001) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x00000002) != 0); }
  bool hasAVX5124VNNIW(void)    { return ((_leaf7sub0edx & 0x00000004) != 0); }
  bool hasAVX5124FMAPS(void)    { return ((_leaf7sub0edx & 0x00000008) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x00000010) != 0); }
  bool hasUINTR(void)           { return ((_leaf7sub0edx & 0x00000020) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x00000040) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x00000080) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x00000100) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x00000200) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x00000400) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x00000800) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x00001000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x00002000) != 0); }
  bool hasSERIALIZE(void)       { return ((_leaf7sub0edx & 0x00004000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x00008000) != 0); }
  bool hasTSXLDTRK(void)        { return ((_leaf7sub0edx & 0x00010000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x00020000) != 0); }
  bool hasPCONFIG(void)         { return ((_leaf7sub0edx & 0x00040000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x00080000) != 0); }
  bool hasIBT(void)             { return ((_leaf7sub0edx & 0x00100000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x00200000) != 0); }
  bool hasAMXBF16(void)         { return ((_leaf7sub0edx & 0x00400000) != 0); }
  bool hasAVX512FP16(void)      { return ((_leaf7sub0edx & 0x00800000) != 0); }
  bool hasAMXTILE(void)         { return ((_leaf7sub0edx & 0x01000000) != 0); }
  bool hasAMXINT8(void)         { return ((_leaf7sub0edx & 0x02000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x04000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x08000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x10000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x20000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x40000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub0edx & 0x80000000) != 0); }

  //  Leaf 7 SubLeaf 1 EAX

  //ol hasRESERVED()            { return ((_leaf7sub1eax & 0x00000001) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1eax & 0x00000002) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1eax & 0x00000004) != 0); }
  bool hasRAOINT(void)          { return ((_leaf7sub1eax & 0x00000008) != 0); }
  bool hasAVXVNNI(void)         { return ((_leaf7sub1eax & 0x00000010) != 0); }
  bool hasAVX512BF16(void)      { return ((_leaf7sub1eax & 0x00000020) != 0); }
  bool hasLASS(void)            { return ((_leaf7sub1eax & 0x00000040) != 0); }
  bool hasCMPCCXADD(void)       { return ((_leaf7sub1eax & 0x00000080) != 0); }
  bool hasARCHPERFMON(void)     { return ((_leaf7sub1eax & 0x00000100) != 0); }   //  Leaf 23h exists
  //ol hasRESERVED()            { return ((_leaf7sub1eax & 0x00000200) != 0); }
  bool hasFAST_MOVSB(void)      { return ((_leaf7sub1eax & 0x00000400) != 0); }
  bool hasFAST_STOSB(void)      { return ((_leaf7sub1eax & 0x00000800) != 0); }
  bool hasFAST_CMPSB(void)      { return ((_leaf7sub1eax & 0x00001000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1eax & 0x00002000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1eax & 0x00004000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1eax & 0x00008000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1eax & 0x00010000) != 0); }
  bool hasFRED(void)            { return ((_leaf7sub1eax & 0x00020000) != 0); }
  bool hasLKGS(void)            { return ((_leaf7sub1eax & 0x00040000) != 0); }
  bool hasWRMSRNS(void)         { return ((_leaf7sub1eax & 0x00080000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1eax & 0x00100000) != 0); }
  bool hasAMXFP16(void)         { return ((_leaf7sub1eax & 0x00200000) != 0); }
  bool hasHRESET(void)          { return ((_leaf7sub1eax & 0x00400000) != 0); }
  bool hasAVXIFMA(void)         { return ((_leaf7sub1eax & 0x00800000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1eax & 0x01000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1eax & 0x02000000) != 0); }
  bool hasLAM(void)             { return ((_leaf7sub1eax & 0x04000000) != 0); }
  bool hasMSRLIST(void)         { return ((_leaf7sub1eax & 0x08000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1eax & 0x10000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1eax & 0x20000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1eax & 0x30000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1eax & 0x40000000) != 0); }

  //  Leaf 7 SubLeaf 1 EDX

  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x00000001) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x00000002) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x00000004) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x00000008) != 0); }
  bool hasAVXVNNIINT8(void)     { return ((_leaf7sub1edx & 0x00000010) != 0); }
  bool hasAVXNECONVERT(void)    { return ((_leaf7sub1edx & 0x00000020) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x00000040) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x00000080) != 0); }
  bool hasAMXCOMPLEX(void)      { return ((_leaf7sub1edx & 0x00000100) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x00000200) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x00000400) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x00000800) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x00001000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x00002000) != 0); }
  bool hasPREFETCHI(void)       { return ((_leaf7sub1edx & 0x00004000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x00008000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x00010000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x00020000) != 0); }
  bool hasCETSSS(void)          { return ((_leaf7sub1edx & 0x00040000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x00080000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x00100000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x00200000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x00300000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x00400000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x01000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x02000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x04000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x08000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x10000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x20000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x40000000) != 0); }
  //ol hasRESERVED()            { return ((_leaf7sub1edx & 0x80000000) != 0); }

  //  Leaf 0xD SubLeaf 1 EAX

  bool hasXSAVEOPT(void)        { return ((_leafDsub0eax & 0x00000001) != 0); }
  bool hasXSAVEC(void)          { return ((_leafDsub0eax & 0x00000002) != 0); }
  bool hasXGETBV_EXC1(void)     { return ((_leafDsub0eax & 0x00000004) != 0); }
  bool hasXSS(void)             { return ((_leafDsub0eax & 0x00000008) != 0); }

  //
  //  AMD features
  //

  //  Leaf 0x80000001 SubLeaf 0 ECX

  bool hasLAHF_LM(void)         { return ((_extf1sub0ecx & 0x00000001) != 0); }   //  dup
  bool hasCMP_LEGACY(void)      { return ((_extf1sub0ecx & 0x00000002) != 0); }   //  dup
  bool hasSVM(void)             { return ((_extf1sub0ecx & 0x00000004) != 0); }   //  dup
  bool hasEXTAPIC(void)         { return ((_extf1sub0ecx & 0x00000008) != 0); }   //  dup
  bool hasCR8_LEGACY(void)      { return ((_extf1sub0ecx & 0x00000010) != 0); }   //  dup
  bool hasABM(void)             { return ((_extf1sub0ecx & 0x00000020) != 0); }   //  dup
  bool hasSSE4A(void)           { return ((_extf1sub0ecx & 0x00000040) != 0); }   //  dup
  bool hasMISALIGNSSE(void)     { return ((_extf1sub0ecx & 0x00000080) != 0); }   //  dup
  bool has3DNOWPREFETCH(void)   { return ((_extf1sub0ecx & 0x00000100) != 0); }   //  dup
  bool hasOSVW(void)            { return ((_extf1sub0ecx & 0x00000200) != 0); }   //  dup
  bool hasIBS(void)             { return ((_extf1sub0ecx & 0x00000400) != 0); }   //  bit 10, reserved
  bool hasXOP(void)             { return ((_extf1sub0ecx & 0x00000800) != 0); }   //  bit 10, reserved
  bool hasSKINIT(void)          { return ((_extf1sub0ecx & 0x00001000) != 0); }   //  dup
  bool hasWDT(void)             { return ((_extf1sub0ecx & 0x00002000) != 0); }   //  dup
  //ol hasRESERVED()            { return ((_extf1sub0ecx & 0x00004000) != 0); }   //  dup
  bool hasLWP(void)             { return ((_extf1sub0ecx & 0x00008000) != 0); }   //  dup
  bool hasFMA4(void)            { return ((_extf1sub0ecx & 0x00010000) != 0); }   //  dup
  bool hasTCE(void)             { return ((_extf1sub0ecx & 0x00020000) != 0); }   //  dup
  //ol hasRESERVED()            { return ((_extf1sub0ecx & 0x00040000) != 0); }   //  bit 18, reserved
  bool hasNODEID_MSR(void)      { return ((_extf1sub0ecx & 0x00080000) != 0); }   //  bit 19, multiprocessor capable
  //ol hasRESERVED()            { return ((_extf1sub0ecx & 0x00100000) != 0); }   //  bit 20, nx bit
  bool hasTBM(void)             { return ((_extf1sub0ecx & 0x00200000) != 0); }   //  bit 21, reserved
  bool hasTOPOEXT(void)         { return ((_extf1sub0ecx & 0x00400000) != 0); }   //  bit 22, extended mmx
  bool hasPERFCTR_CORE(void)    { return ((_extf1sub0ecx & 0x00800000) != 0); }   //  dup
  bool hasPERFCTR_NB(void)      { return ((_extf1sub0ecx & 0x01000000) != 0); }   //  dup
  //ol hasRESERVED()            { return ((_extf1sub0ecx & 0x02000000) != 0); }   //  bit 25, fxsave/fxrstor optimizations 
  bool hasDBX(void)             { return ((_extf1sub0ecx & 0x04000000) != 0); }   //  bit 26, gigabyte pages
  bool hasPERFTSC(void)         { return ((_extf1sub0ecx & 0x08000000) != 0); }   //  bit 27, rdtscp instruction
  bool hasPCX_L2I(void)         { return ((_extf1sub0ecx & 0x10000000) != 0); }   //  bit 28, reserved
  bool hasMONITOR(void)         { return ((_extf1sub0ecx & 0x20000000) != 0); }   //  bit 29, long mode
  bool hasADDR_MASK_EXT(void)   { return ((_extf1sub0ecx & 0x40000000) != 0); }   //  bit 30, extended 3DNow!
  //ol hasRESERVED()            { return ((_extf1sub0ecx & 0x80000000) != 0); }   //  bit 31, 3DNow!

  //  Leaf 0x80000001 SubLeaf 0 EDX
  //
  //  This returns extended feature flags in EDX and ECX. Bits 0 through 9,
  //  12 through 17, 23, and 24 of EDX are duplicates of EDX from the EAX=1
  //  leaf.

  //ol hasDUPLICATE(void)       { return ((_extf1sub0edx & 0x00000001) != 0); }   //  dup
  //ol hasDUPLICATE(void)       { return ((_extf1sub0edx & 0x00000002) != 0); }   //  dup
  //ol hasDUPLICATE(void)       { return ((_extf1sub0edx & 0x00000004) != 0); }   //  dup
  //ol hasDUPLICATE(void)       { return ((_extf1sub0edx & 0x00000008) != 0); }   //  dup
  //ol hasDUPLICATE(void)       { return ((_extf1sub0edx & 0x00000010) != 0); }   //  dup
  //ol hasDUPLICATE(void)       { return ((_extf1sub0edx & 0x00000020) != 0); }   //  dup
  //ol hasDUPLICATE(void)       { return ((_extf1sub0edx & 0x00000040) != 0); }   //  dup
  //ol hasDUPLICATE(void)       { return ((_extf1sub0edx & 0x00000080) != 0); }   //  dup
  //ol hasDUPLICATE(void)       { return ((_extf1sub0edx & 0x00000100) != 0); }   //  dup
  //ol hasDUPLICATE(void)       { return ((_extf1sub0edx & 0x00000200) != 0); }   //  dup
  //ol hasRESERVED()            { return ((_extf1sub0edx & 0x00000400) != 0); }   //  bit 10
  bool hasSYSCALL(void)         { return ((_extf1sub0edx & 0x00000800) != 0); }   //  bit 11
  //ol hasDUPLICATE(void)       { return ((_extf1sub0edx & 0x00001000) != 0); }   //  dup
  //ol hasDUPLICATE(void)       { return ((_extf1sub0edx & 0x00002000) != 0); }   //  dup
  //ol hasDUPLICATE(void)       { return ((_extf1sub0edx & 0x00004000) != 0); }   //  dup
  //ol hasDUPLICATE(void)       { return ((_extf1sub0edx & 0x00008000) != 0); }   //  dup
  //ol hasDUPLICATE(void)       { return ((_extf1sub0edx & 0x00010000) != 0); }   //  dup
  //ol hasDUPLICATE(void)       { return ((_extf1sub0edx & 0x00020000) != 0); }   //  dup
  //ol hasRESERVED()            { return ((_extf1sub0edx & 0x00040000) != 0); }   //  bit 18, reserved
  bool hasMP(void)              { return ((_extf1sub0edx & 0x00080000) != 0); }   //  bit 19, multiprocessor capable
  bool hasNX(void)              { return ((_extf1sub0edx & 0x00100000) != 0); }   //  bit 20, nx bit
  //ol hasRESERVED()            { return ((_extf1sub0edx & 0x00200000) != 0); }   //  bit 21, reserved
  bool hasMMXEXT(void)          { return ((_extf1sub0edx & 0x00400000) != 0); }   //  bit 22, extended mmx
  //ol hasDUPLICATE(void)       { return ((_extf1sub0edx & 0x00800000) != 0); }   //  dup
  //ol hasDUPLICATE(void)       { return ((_extf1sub0edx & 0x01000000) != 0); }   //  dup
  bool hasFXSR_OPT(void)        { return ((_extf1sub0edx & 0x02000000) != 0); }   //  bit 25, fxsave/fxrstor optimizations 
  bool hasPDPE1GB(void)         { return ((_extf1sub0edx & 0x04000000) != 0); }   //  bit 26, gigabyte pages
  bool hasRDTSCPM(void)         { return ((_extf1sub0edx & 0x08000000) != 0); }   //  bit 27, rdtscp instruction
  //ol hasRESERVED()            { return ((_extf1sub0edx & 0x10000000) != 0); }   //  bit 28, reserved
  bool hasLM(void)              { return ((_extf1sub0edx & 0x20000000) != 0); }   //  bit 29, long mode
  bool has3DNOWEXT(void)        { return ((_extf1sub0edx & 0x40000000) != 0); }   //  bit 30, extended 3DNow!
  bool has3DNOW(void)           { return ((_extf1sub0edx & 0x80000000) != 0); }   //  bit 31, 3DNow!

private:
  uint32  _maxLeaf    = 0;
  uint32  _maxExtLeaf = 0;

  bool    isLeafValid(uint32 leaf)          { return               leaf  <= _maxLeaf;    }
  bool    isExtendedLeafValid(uint32 leaf)  { return (0x80000000 | leaf) <= _maxExtLeaf; }

  void    loadProcessorFlags(void);

  void    decodeProcessorOrigin(void);
  void    decodeProcessorName(void);

  void    decodeProcessorModelVendorStepID(void);

  bool    decodeIntelCacheTopology(void);
  bool    decodeAMDCacheTopology(void);
  bool    decodeAMDLegacyCacheTopology(void);

public:
  char    _processorName[48]   = {0};
  char    _processorOrigin[16] = {0};

  uint32  _modelID  = 0;
  uint32  _familyID = 0;
  uint32  _stepping = 0;

  bool    _isIntel = false;
  bool    _isAMD = false;

private:
  uint64  _xcr0;

  uint32  _leaf1sub0eax, _leaf1sub0ebx, _leaf1sub0ecx, _leaf1sub0edx;

  uint32  _leaf7sub0eax, _leaf7sub0ebx, _leaf7sub0ecx, _leaf7sub0edx;
  uint32  _leaf7sub1eax, _leaf7sub1ebx, _leaf7sub1ecx, _leaf7sub1edx;
  uint32  _leafDsub0eax, _leafDsub0ebx, _leafDsub0ecx, _leafDsub0edx;

  uint32  _extf1sub0eax, _extf1sub0ebx, _extf1sub0ecx, _extf1sub0edx;
};

}

#endif  //  MERYLUTIL_SYSTEM_CPUID_H
