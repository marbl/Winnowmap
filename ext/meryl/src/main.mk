MODULE       :=    meryl
TARGET       := libmeryl.a
SOURCES      := utility/src/utility/edlib.C \
                \
                utility/src/utility/files.C \
                utility/src/utility/files-buffered.C \
                utility/src/utility/files-compressed.C \
                utility/src/utility/files-fasta-fastq.C \
                utility/src/utility/files-memoryMapped.C \
                \
                utility/src/utility/logging.C \
                \
                utility/src/utility/strings.C \
                \
                utility/src/utility/system.C \
                utility/src/utility/system-stackTrace.C \
                \
                utility/src/utility/sequence.C \
                \
                utility/src/utility/types.C \
                \
                utility/src/utility/kmers-exact.C \
                utility/src/utility/kmers-files.C \
                utility/src/utility/kmers-histogram.C \
                utility/src/utility/kmers-reader.C \
                utility/src/utility/kmers-writer-block.C \
                utility/src/utility/kmers-writer-stream.C \
                utility/src/utility/kmers-writer.C \
                utility/src/utility/kmers.C \
                \
                utility/src/utility/bits.C \
                utility/src/utility/bits-wordArray.C \
                \
                utility/src/utility/hexDump.C \
                utility/src/utility/md5.C \
                utility/src/utility/mt19937ar.C \
                utility/src/utility/speedCounter.C \
                utility/src/utility/sweatShop.C \
                \
                utility/src/utility/runtime.C

ifeq (${BUILDSTACKTRACE}, 1)
SOURCES      += utility/src/utility/libbacktrace/atomic.c \
                utility/src/utility/libbacktrace/backtrace.c \
                utility/src/utility/libbacktrace/dwarf.c \
                utility/src/utility/libbacktrace/elf.c \
                utility/src/utility/libbacktrace/fileline.c \
                utility/src/utility/libbacktrace/mmap.c \
                utility/src/utility/libbacktrace/mmapio.c \
                utility/src/utility/libbacktrace/posix.c \
                utility/src/utility/libbacktrace/print.c \
                utility/src/utility/libbacktrace/simple.c \
                utility/src/utility/libbacktrace/sort.c \
                utility/src/utility/libbacktrace/state.c \
                utility/src/utility/libbacktrace/unknown.c
endif

SRC_INCDIRS  := . \
                meryl \
                utility

SUBMAKEFILES := meryl/meryl.mk \
                meryl-analyze/meryl-analyze.mk \
                meryl-simple/meryl-simple.mk \
                meryl-import/meryl-import.mk \
                meryl-lookup/meryl-lookup.mk

ifeq ($(BUILDTESTS), 1)
SUBMAKEFILES += tests/merylCountArrayTest.mk \
                tests/merylExactLookupTest.mk
endif
