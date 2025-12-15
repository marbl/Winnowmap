
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
#include "htslib/hts/sam.h"

#include <assert.h>

//using namespace merylutil;


int
main(int argc, char **argv) {


    const char *qname = "q1";
    const uint32_t cigar[] = { 6 << BAM_CIGAR_SHIFT | BAM_CMATCH, 2 << BAM_CIGAR_SHIFT | BAM_CINS, 2 << BAM_CIGAR_SHIFT | BAM_CMATCH };
    const char *seq = "TGGACTACGA";
    const char *qual = "DBBBB+=7=0";
    const char *temp_fname = "testHTSbam.bam";

    int r;

    //htsFile *writer, *reader;
    //bam1_t *w_bam, *r_bam;
    kstring_t ks = KS_INITIALIZE;




    // open file for writing
    htsFile *writer = hts_open(temp_fname, "wb");
    assert(writer != NULL);

    // write header then align
    if (1) {
    sam_hdr_t *w_header = bam_hdr_init();
    assert(w_header != NULL);

    r = sam_hdr_add_line(w_header, "SQ", "SN", "t1", "LN", "5001", NULL);
    assert(r == 0);

    r = sam_hdr_write(writer, w_header);
    assert(r == 0);

    // write alignments
    bam1_t *w_bam = bam_init1();
    assert(w_bam != NULL);

    r = bam_set1(w_bam, strlen(qname), qname,
                 BAM_FREVERSE, 0, 1000, 42,
                 sizeof(cigar) / 4, cigar, 0, 2000, 3000,
                 strlen(seq), seq, qual, 64);
    assert(r >= 0);

    r = sam_write1(writer, w_header, w_bam);
    assert(r >= 0);

    r = sam_write1(writer, w_header, w_bam);
    assert(r >= 0);

    bam_destroy1(w_bam);
    sam_hdr_destroy(w_header);
    }



    if (0) {
    sam_hdr_t *w_header = bam_hdr_init();
    assert(w_header != NULL);

    r = sam_hdr_add_line(w_header, "SQ", "SN", "t1", "LN", "5001", NULL);
    assert(r == 0);

    r = sam_hdr_write(writer, w_header);
    assert(r == 0);

    // write alignments
    bam1_t *w_bam = bam_init1();
    assert(w_bam != NULL);

    r = bam_set1(w_bam, strlen(qname), qname,
                 BAM_FREVERSE, 0, 1000, 42,
                 sizeof(cigar) / 4, cigar, 0, 2000, 3000,
                 strlen(seq), seq, qual, 64);
    assert(r >= 0);

    r = sam_write1(writer, w_header, w_bam);
    assert(r >= 0);

    bam_destroy1(w_bam);
    sam_hdr_destroy(w_header);
    }





    // close file
    r = hts_close(writer);
    assert(r == 0);



#if 0
    // open file for reading
    htsFile *reader = hts_open(temp_fname, "rb");
    assert(reader != NULL);

    // read header
    sam_hdr_t *r_header = sam_hdr_read(reader);
    assert(r_header != NULL);

    r = sam_hdr_find_tag_id(r_header, "SQ", NULL, NULL, "SN", &ks);
    assert(r == 0);
    assert(strcmp(ks_c_str(&ks), "t1") == 0);
    assert(r_header->n_targets == 1);
    assert(strcmp(r_header->target_name[0], "t1") == 0);
    assert(r_header->target_len[0] == 5000);

    // read alignments
    bam1_t *r_bam = bam_init1();
    assert(r_bam != NULL);

    r = sam_read1(reader, r_header, r_bam);
    assert(r >= 0);
    assert(strcmp(bam_get_qname(r_bam), qname) == 0);
    assert(r_bam->core.n_cigar == sizeof(cigar) / 4);
    assert(memcmp(bam_get_cigar(r_bam), cigar, sizeof(cigar)) == 0);
    assert(r_bam->core.l_qseq == strlen(seq));

    r = sam_read1(reader, r_header, r_bam);
    assert(r < 0);

    bam_destroy1(r_bam);

    // close file
    r = hts_close(reader);
    assert(r == 0);

    sam_hdr_destroy(r_header);
#endif

    ks_free(&ks);

  return 0;
}

