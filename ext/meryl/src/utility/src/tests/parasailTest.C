

#include "system.H"
#include "types.H"
#include "files.H"
#include "sequence.H"

#include "align-parasail-driver.H"

#include "parasail.h"
#include "parasail/cpuid.h"



void
checkEncodeDecodeBase(void) {
  //fprintf(stderr, "A -> %2u -> %c\n", encode2bitBase('A'), decode2bitBase(0));
  assert(encode2bitBase('A') == 0);
  assert(decode2bitBase(0) == 'A');

  //fprintf(stderr, "C -> %2u -> %c\n", encode2bitBase('C'), decode2bitBase(1));
  assert(encode2bitBase('C') == 1);
  assert(decode2bitBase(1) == 'C');

  //fprintf(stderr, "G -> %2u -> %c\n", encode2bitBase('G'), decode2bitBase(3));
  assert(encode2bitBase('G') == 3);
  assert(decode2bitBase(3) == 'G');

  //fprintf(stderr, "T -> %2u -> %c\n", encode2bitBase('T'), decode2bitBase(2));
  assert(encode2bitBase('T') == 2);
  assert(decode2bitBase(2) == 'T');
};



void
checkVectorSupport(void) {
  fprintf(stderr, "avx512vbmi:  %d\n", parasail_can_use_avx512vbmi());
  fprintf(stderr, "avx512bw:    %d\n", parasail_can_use_avx512bw());
  fprintf(stderr, "avx512f:     %d\n", parasail_can_use_avx512f());
  fprintf(stderr, "avx2:        %d\n", parasail_can_use_avx2());
  fprintf(stderr, "sse41:       %d\n", parasail_can_use_sse41());
  fprintf(stderr, "sse2:        %d\n", parasail_can_use_sse2());
}


void
printCigar(parasail_result_t *result,
           parasail_matrix_t *matrix,
           char const *seqA, uint32 lenA,
           char const *seqB, uint32 lenB) {
  parasail_cigar_t  *cigar = parasail_result_get_cigar(result, seqA, lenA, seqB, lenB, matrix);

  fprintf(stderr, "A: %6d-%6d\n", cigar->beg_query, parasail_result_get_end_query(result));
  fprintf(stderr, "B: %6d-%6d\n", cigar->beg_ref,   parasail_result_get_end_ref  (result));

  if (cigar->len <= 30) {
    fprintf(stderr, "CIGAR:");
    for (uint32 cc=0; cc<cigar->len; cc++)
      fprintf(stderr, " %d%c", parasail_cigar_decode_len(cigar->seq[cc]), parasail_cigar_decode_op(cigar->seq[cc]));
    fprintf(stderr, "\n");
  }

  else {
    fprintf(stderr, "CIGAR: %d elements.\n", cigar->len);
  }

  parasail_cigar_free(cigar);
}


void
explicitCallParasail(char const *seqA, uint32 lenA,
                     char const *seqB, uint32 lenB) {

  //  

  parasail_matrix_t  *matrix = parasail_matrix_create("ACGTN", 2, -1);
  parasail_result_t  *result = nullptr;
  parasail_cigar_t   *cigar  = nullptr;
  double              bgn    = 0.0;

  //  Gap 'open' and 'extend' are passed as positive penalty values.
  //  'open' must be at least as big as 'extend'.
  //  The first sequence is the 'query', the second is the 'database'.
  //
  //  sg       -- free gaps at all ends
  //  sg_qb_de -- free gaps at the start of s1 and the end of s2
  //  sg_qe_db -- free gaps at the end of s1 and the start of s2

  fprintf(stderr, "\n--\n");
  fprintf(stderr, "non-vector A-vs-B\n");
  bgn = getTime();
  result = parasail_sg_qb_de_trace(seqA, lenA, seqB, lenB, 2, 1, matrix);
  fprintf(stderr, "%.3f seconds  score %d\n", getTime() - bgn, parasail_result_get_score(result));
  printCigar(result, matrix, seqA, lenA, seqB, lenB);
  parasail_result_free(result);

  fprintf(stderr, "\n--\n");
  fprintf(stderr, "non-vector B-vs-A\n");
  bgn = getTime();
  result = parasail_sg_qb_de_trace(seqB, lenB, seqA, lenA, 2, 1, matrix);
  fprintf(stderr, "%.3f seconds  score %d\n", getTime() - bgn, parasail_result_get_score(result));
  printCigar(result, matrix, seqA, lenA, seqB, lenB);
  parasail_result_free(result);

#if 0
  fprintf(stderr, "\n--\n");
  fprintf(stderr, "vector\n");
  bgn = getTime();
  result = parasail_sg_qb_de_trace_striped_16(seqA, lenA, seqB, lenB, 2, 1, matrix);
  fprintf(stderr, "%.3f seconds  score %d\n", getTime() - bgn, parasail_result_get_score(result));
  printCigar(result, matrix, seqA, lenA, seqB, lenB);
  parasail_result_free(result);

  fprintf(stderr, "\n--\n");
  fprintf(stderr, "avx2\n");
  bgn = getTime();
  result = parasail_sg_qb_de_trace_striped_avx2_256_16(seqA, lenA, seqB, lenB, 2, 1, matrix);
  fprintf(stderr, "%.3f seconds  score %d\n", getTime() - bgn, parasail_result_get_score(result));
  printCigar(result, matrix, seqA, lenA, seqB, lenB);
  parasail_result_free(result);

  fprintf(stderr, "\n--\n");
  fprintf(stderr, "sse2\n");
  bgn = getTime();
  result = parasail_sg_qb_de_trace_striped_sse2_128_16(seqA, lenA, seqB, lenB, 2, 1, matrix);
  fprintf(stderr, "%.3f seconds  score %d\n", getTime() - bgn, parasail_result_get_score(result));
  printCigar(result, matrix, seqA, lenA, seqB, lenB);
  parasail_result_free(result);

  fprintf(stderr, "\n--\n");
  fprintf(stderr, "sse4.1\n");
  bgn = getTime();
  result = parasail_sg_qb_de_trace_striped_sse41_128_16(seqA, lenA, seqB, lenB, 2, 1, matrix);
  fprintf(stderr, "%.3f seconds  score %d\n", getTime() - bgn, parasail_result_get_score(result));
  printCigar(result, matrix, seqA, lenA, seqB, lenB);
  parasail_result_free(result);
#endif

  parasail_matrix_free(matrix);
}




int
main(int argc, char **argv) {
  char const *seqA = nullptr;   uint32  lenA = 0;   dnaSeq dseqA;
  char const *seqB = nullptr;   uint32  lenB = 0;   dnaSeq dseqB;

  bool  alignDovetail  = false;
  bool  alignContained = false;
  bool  explicitTests  = false;

  checkEncodeDecodeBase();
  checkVectorSupport();

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if      (strcmp(argv[arg], "-a") == 0) {
      seqA = argv[++arg];
      lenA = strlen(seqA);
    }

    else if (strcmp(argv[arg], "-b") == 0) {
      seqB = argv[++arg];
      lenB = strlen(seqB);
    }

    else if (strcmp(argv[arg], "-A") == 0) {
      dnaSeqFile  fileA(argv[++arg]);

      fileA.loadSequence(dseqA);

      seqA = dseqA.bases();
      lenA = strlen(seqA);
    }

    else if (strcmp(argv[arg], "-B") == 0) {
      dnaSeqFile  fileB(argv[++arg]);

      fileB.loadSequence(dseqB);

      seqB = dseqB.bases();
      lenB = strlen(seqB);
    }

    else if (strcmp(argv[arg], "-C") == 0) {
      alignContained = true;
    }

    else if (strcmp(argv[arg], "-D") == 0) {
      alignDovetail = true;
    }

    else if (strcmp(argv[arg], "-E") == 0) {
      explicitTests = true;
    }

    else {
      err++;
    }

    arg++;
  }

  if ((alignDovetail == false) && (alignContained == false))  err++;

  if ((seqA == nullptr) || (lenA == 0))   err++;
  if ((seqB == nullptr) || (lenB == 0))   err++;

  if (err > 0) {
    fprintf(stderr, "usage: %s [-C | -D] -a <bases> -b <bases>\n", argv[0]);
    fprintf(stderr, "       %*s [-C | -D] -A <file>  -B <file>\n", (int)strlen(argv[0]), "");
    fprintf(stderr, "       %*s [-E]\n", (int)strlen(argv[0]), "");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -C  align A as the container of B:\n");
    fprintf(stderr, "        A -------------\n");
    fprintf(stderr, "        B     -----\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -D  align A with a dovetail overlap to B\n");
    fprintf(stderr, "        A -------------\n");
    fprintf(stderr, "        B     -------------\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -E  run tests calling Parasail directly; assumes\n");
    fprintf(stderr, "      sequences are dovetail as in -D.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "Testing A: '%20.20s....' length %u against\n", seqA, lenA);
  fprintf(stderr, "        B: '%20.20s....' length %u\n",         seqB, lenB);
  fprintf(stderr, "\n");

  if (explicitTests)
    explicitCallParasail(seqA, lenA, seqB, lenB);

  {
    parasailLib  pl(2, -1, 2, 1);

    fprintf(stderr, "\n--\n");
    fprintf(stderr, "parasailLib\n");

    if (alignDovetail)
      pl.alignDovetail(seqA, lenA, 0, lenA,
                       seqB, lenB, 0, lenB, true);

    if (alignContained)
      pl.alignContained(seqA, lenA, 0, lenA,
                        seqB, lenB, 0, lenB, true);

    fprintf(stderr, "----------\n");
    fprintf(stderr, "A: %5u-%5u\n",         pl.bgnA(), pl.endA());
    fprintf(stderr, "B: %5u-%5u  %.4f%%\n", pl.bgnB(), pl.endB(), pl.percentIdentity());
  }

  return(0);
}
