MODULE       := meryl
VERSION      := snapshot 1.3
VERSION_H    := utility/src/version.H

TARGET       := libmeryl.a
SOURCES      := utility/src/align/align-ksw2-driver.C \
                utility/src/align/align-ksw2-extz.C \
                utility/src/align/align-ksw2-extz2-sse.C \
                utility/src/align/align-parasail-driver.C \
                utility/src/align/align-ssw-driver.C \
                utility/src/align/align-ssw.C \
                utility/src/align/edlib.C \
                \
                utility/src/bits/fibonacci-v1.C \
                utility/src/bits/hexDump-v1.C \
                utility/src/bits/stuffedBits-v1-binary.C \
                utility/src/bits/stuffedBits-v1-bits.C \
                utility/src/bits/stuffedBits-v1-delta.C \
                utility/src/bits/stuffedBits-v1-gamma.C \
                utility/src/bits/stuffedBits-v1-golomb.C \
                utility/src/bits/stuffedBits-v1-omega.C \
                utility/src/bits/stuffedBits-v1-unary.C \
                utility/src/bits/stuffedBits-v1-zeckendorf.C \
                utility/src/bits/stuffedBits-v1.C \
                utility/src/bits/wordArray-v1.C \
                \
                utility/src/datastructures/keyAndValue-v1.C \
                \
                utility/src/datastructures/regex-v1.C \
                utility/src/datastructures/regex-v2.C \
                utility/src/datastructures/regex-v2-build.C \
                utility/src/datastructures/regex-v2-convert.C \
                utility/src/datastructures/regex-v2-match.C \
                utility/src/datastructures/regex-v2-parse.C \
                \
                utility/src/datastructures/splitToWords-v1.C \
                utility/src/datastructures/stringList-v1.C \
                utility/src/datastructures/strings-v1.C \
                utility/src/datastructures/types-v1.C \
                \
                utility/src/files/accessing-v1.C \
                utility/src/files/buffered-v1-reading.C \
                utility/src/files/buffered-v1-writing.C \
                utility/src/files/compressed-v1-reading.C \
                utility/src/files/compressed-v1-writing.C \
                utility/src/files/compressed-v1.C \
                utility/src/files/fasta-fastq-v1.C \
                utility/src/files/files-v1.C \
                utility/src/files/memoryMapped-v1.C \
                utility/src/files/readLine-v0.C \
                utility/src/files/readLine-v1.C \
                utility/src/files/reading-v1.C \
                utility/src/files/writing-v1.C \
                \
                utility/src/kmers-v1/kmers-exact.C \
                utility/src/kmers-v1/kmers-files.C \
                utility/src/kmers-v1/kmers-histogram.C \
                utility/src/kmers-v1/kmers-histogram-ploidy.C \
                utility/src/kmers-v1/kmers-reader.C \
                utility/src/kmers-v1/kmers-writer-block.C \
                utility/src/kmers-v1/kmers-writer-stream.C \
                utility/src/kmers-v1/kmers-writer.C \
                utility/src/kmers-v1/kmers.C \
                \
                utility/src/math/md5-v1.C \
                utility/src/math/mt19937ar-v1.C \
                utility/src/math/sampledDistribution-v1.C \
                \
                utility/src/sequence/dnaSeq-v1.C \
                utility/src/sequence/bufSeqFile-v1.C \
                utility/src/sequence/htsSeqFile-v1.C \
                utility/src/sequence/sequence-v1.C \
                \
                utility/src/system/logging-v1.C \
                utility/src/system/runtime-v1.C \
                utility/src/system/speedCounter-v1.C \
                utility/src/system/sweatShop-v1.C \
                utility/src/system/system-stackTrace-v1.C \
                utility/src/system/system-v1.C \
                utility/src/system/time-v1.C

#SOURCES      += utility/src/kmers-v2/kmers-exact.C \
#                utility/src/kmers-v2/kmers-files.C \
#                utility/src/kmers-v2/kmers-histogram.C \
#                utility/src/kmers-v2/kmers-reader-dump.C \
#                utility/src/kmers-v2/kmers-reader.C \
#                utility/src/kmers-v2/kmers-writer-block.C \
#                utility/src/kmers-v2/kmers-writer-stream.C \
#                utility/src/kmers-v2/kmers-writer.C \
#                utility/src/kmers-v2/kmers.C

SOURCES      += utility/src/htslib/hts/bcf_sr_sort.c \
                utility/src/htslib/hts/bgzf.c \
                utility/src/htslib/hts/errmod.c \
                utility/src/htslib/hts/faidx.c \
                utility/src/htslib/hts/header.c \
                utility/src/htslib/hts/hfile.c \
                utility/src/htslib/hts/hfile_libcurl.c \
                utility/src/htslib/hts/hfile_s3.c \
                utility/src/htslib/hts/hts.c \
                utility/src/htslib/hts/hts_expr.c \
                utility/src/htslib/hts/hts_os.c \
                utility/src/htslib/hts/kfunc.c \
                utility/src/htslib/hts/kstring.c \
                utility/src/htslib/hts/md5.c \
                utility/src/htslib/hts/multipart.c \
                utility/src/htslib/hts/probaln.c \
                utility/src/htslib/hts/realn.c \
                utility/src/htslib/hts/regidx.c \
                utility/src/htslib/hts/region.c \
                utility/src/htslib/hts/sam.c \
                utility/src/htslib/hts/synced_bcf_reader.c \
                utility/src/htslib/hts/tbx.c \
                utility/src/htslib/hts/textutils.c \
                utility/src/htslib/hts/thread_pool.c \
                utility/src/htslib/hts/vcf.c \
                utility/src/htslib/hts/vcf_sweep.c \
                utility/src/htslib/hts/vcfutils.c \
                utility/src/htslib/cram/cram_codecs.c \
                utility/src/htslib/cram/cram_decode.c \
                utility/src/htslib/cram/cram_encode.c \
                utility/src/htslib/cram/cram_external.c \
                utility/src/htslib/cram/cram_index.c \
                utility/src/htslib/cram/cram_io.c \
                utility/src/htslib/cram/cram_stats.c \
                utility/src/htslib/cram/mFILE.c \
                utility/src/htslib/cram/open_trace_file.c \
                utility/src/htslib/cram/pooled_alloc.c \
                utility/src/htslib/cram/string_alloc.c \
                utility/src/htslib/htscodecs/arith_dynamic.c \
                utility/src/htslib/htscodecs/fqzcomp_qual.c \
                utility/src/htslib/htscodecs/htscodecs.c \
                utility/src/htslib/htscodecs/pack.c \
                utility/src/htslib/htscodecs/rANS_static.c \
                utility/src/htslib/htscodecs/rANS_static32x16pr.c \
                utility/src/htslib/htscodecs/rANS_static32x16pr_neon.c \
                utility/src/htslib/htscodecs/rANS_static4x16pr.c \
                utility/src/htslib/htscodecs/rle.c \
                utility/src/htslib/htscodecs/tokenise_name3.c \
                utility/src/htslib/htscodecs/utils.c


ifeq (${BUILDSTACKTRACE}, 1)
SOURCES      += utility/src/system/libbacktrace/atomic.c \
                utility/src/system/libbacktrace/backtrace.c \
                utility/src/system/libbacktrace/dwarf.c \
                utility/src/system/libbacktrace/elf.c \
                utility/src/system/libbacktrace/fileline.c \
                utility/src/system/libbacktrace/mmap.c \
                utility/src/system/libbacktrace/mmapio.c \
                utility/src/system/libbacktrace/posix.c \
                utility/src/system/libbacktrace/print.c \
                utility/src/system/libbacktrace/simple.c \
                utility/src/system/libbacktrace/sort.c \
                utility/src/system/libbacktrace/state.c \
                utility/src/system/libbacktrace/unknown.c
endif

SRC_INCDIRS  := . \
                meryl \
                utility/src

SYS_INCDIRS  += $(shell pkg-config --cflags-only-I openssl libcurl liblzma | sed s:-I/:/:g)
LDFLAGS      += $(shell pkg-config --libs-only-L   openssl libcurl liblzma)
LDLIBS       += $(shell pkg-config --libs-only-l   openssl libcurl liblzma) -lz -lbz2

SUBMAKEFILES := meryl/meryl.mk \
                meryl-analyze/meryl-analyze.mk \
                meryl-simple/meryl-simple.mk \
                meryl-import/meryl-import.mk \
                meryl-lookup/meryl-lookup.mk \
                meryl-lookup/position-lookup.mk

#SUBMAKEFILES += meryl2/meryl.mk \
#                meryl2-analyze/meryl-analyze.mk \
#                meryl2-ploidy/meryl-ploidy.mk \
#                meryl2-simple/meryl-simple.mk \
#                meryl2-import/meryl-import.mk \
#                meryl2-lookup/meryl-lookup.mk

ifeq ($(BUILDTESTS), 1)
SUBMAKEFILES += tests/merylCountArrayTest.mk \
                tests/merylExactLookupTest.mk \
                tests/matchTokenTest.mk
endif
