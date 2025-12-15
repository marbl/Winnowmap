//
//  mcaSize       = sizeof(merylCountArray)  == 80
//  ptrSize       = sizeof(uint64 *)         == 8
//  segSize       = (as above)               ~= 64 KB
//
//  mersPerPrefix = nKmers / nPrefix+1
//  mersPerSeg    = 8 * segSize / 2*merSize-prefixSize
//
//  nSegPerPrefix = mersPerPrefix / mersPerSeg
//
//             basic structure             pointers to data                         data
//           |-----------------|   |------------------------------|   |------------------------------|
//           |                 |   |                              |   |                              |
//  memory = (mcaSize * nPrefix) + (ptrSize * nPrefix * nSegPrefix) + (segSize * nPrefix * nSegPrefix)
//
//  Solving for nKmers
//
//  memory = (mcaSize * nPrefix) + nSegPerPrefix * (ptrSize * nPrefix + segSize * nPrefix)
//
//  memory - (mcaSize * nPrefix) = nSegPerPrefix * (ptrSize * nPrefix + segSize * nPrefix)
//
//  nSegPerPrefix = (memory - mcaSize * nPrefix) / (ptrSize * nPrefix + nPrefix * segSize)
//
//  (nKmers / nPrefix+1) / mersPerSeg = (memory - mcaSize * nPrefix) / (ptrSize * nPrefix + segSize * nPrefix)
//   nKmers / nPrefix+1  = mersPerSeg * (memory - mcaSize * nPrefix) / (ptrSize * nPrefix + segSize * nPrefix)
#if 0
uint64
findMaxInputSizeForMemorySize(uint32 merSize, uint64 memSize) {
  uint64  mcaSize = sizeof(merylCountArray);
  uint64  ptrSize = sizeof(uint64 *);
  uint64  segSize = SEGMENT_SIZE * 1024;

  //  Free variable - prefixSize (wp)

  fprintf(stderr, "For memory size %lu bytes\n", memSize);
  fprintf(stderr, "                                                  |---------memory-breakdown--------\n");
  fprintf(stderr, "pBits    nKmers     nPrefix   nSegment   mers/seg  structure   pointers         data\n");
  fprintf(stderr, "----- ---------- ---------- ---------- ---------- ---------- ---------- ------------\n");

  for (uint32 wp=1; wp < 2*merSize; wp++) {
    uint64  nPrefix       = (uint64)1 << wp;

    if (mcaSize * nPrefix > memSize)
      break;

    //  Compute how many mer suffixes we can fit in a segment of memory.

    uint64  mersPerSeg = (segSize * 8) / (2 * merSize - wp);

    //  Compute the number of segments we can fit in memeory.  Each prefix must have a segment.  We
    //  don't actually know how many pointers we'll need -- it depends on the number of prefixes and
    //  number of segments -- so we assume there are 64 per prefix.
    //
    uint64  nSeg = (memSize - mcaSize * nPrefix - ptrSize * 64 * nPrefix) / segSize;

    if (nSeg < nPrefix)
      nSeg = nPrefix;

    //  Assuming some segment loading average, compute the number of kmers we can fit.  Each prefix
    //  has a 2/3 full segment, then any left over segments are 100% full.

    uint64  nKmers = nPrefix * mersPerSeg * 0.666;

    if (nPrefix < nSeg)
      nKmers += (nSeg - nPrefix) * mersPerSeg;

    //  For flavoring, show the memory breakdown.  nSegPrefix underflows, so needs to be included directly.

    uint64  ptrPerPrefix = 64; //32 * (nSeg / nPrefix + 1);  //  Somewhat incorrect.  The basic allocation is 32 pointers, then it doubles.

    uint64  basicMem     = mcaSize * nPrefix;                                    //  Size of the merCountArray structure itself.
    uint64  ptrMem       = ptrSize * nPrefix * ptrPerPrefix;                     //  Additional allocation for pointers to segments.
    uint64  dataMem      = segSize * nSeg;                                       //  Additional allocation of segments.

    if (basicMem + ptrMem + dataMem > 4 * memSize)
      break;

    fprintf(stderr, "%5u %10lu %10lu %10lu %10lu %10lu %10lu %12lu%s\n",
            wp, nKmers,
            nPrefix, nSeg,
            mersPerSeg,
            basicMem, ptrMem, dataMem,
            (basicMem + ptrMem + dataMem < memSize) ? "" : " INVALID");
  }

  exit(0);
}
#endif

