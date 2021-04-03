
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.

ifeq "$(strip ${DESTDIR})" ""
  DESTDIR      :=
endif

ifeq "$(strip ${PREFIX})" ""
  ifeq "$(strip ${DESTDIR})" ""
    PREFIX     := $(realpath ..)
  else
    PREFIX     := /canu
  endif
endif

ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := $(DESTDIR)$(PREFIX)/$(OSTYPE)-$(MACHINETYPE)/obj
endif

ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := $(DESTDIR)$(PREFIX)/$(OSTYPE)-$(MACHINETYPE)
endif

TARGET       := libcanu.a

SOURCES      := utility/runtime.C \
                \
                utility/edlib.C \
                \
                utility/files.C \
                utility/files-buffered.C \
                utility/files-compressed.C \
                utility/files-memoryMapped.C \
                \
                utility/logging.C \
                \
                utility/strings.C \
                \
                utility/system.C \
                utility/system-stackTrace.C \
                \
                utility/sequence.C \
                \
                utility/types.C \
                \
                utility/kmers-exact.C \
                utility/kmers-files.C \
                utility/kmers-histogram.C \
                utility/kmers-reader.C \
                utility/kmers-writer-block.C \
                utility/kmers-writer-stream.C \
                utility/kmers-writer.C \
                utility/kmers.C \
                \
                utility/bits.C \
                \
                utility/hexDump.C \
                utility/md5.C \
                utility/mt19937ar.C \
                utility/objectStore.C \
                utility/speedCounter.C \
                utility/sweatShop.C


ifeq (${BUILDSTACKTRACE}, 1)
SOURCES      += utility/libbacktrace/atomic.c \
                utility/libbacktrace/backtrace.c \
                utility/libbacktrace/dwarf.c \
                utility/libbacktrace/elf.c \
                utility/libbacktrace/fileline.c \
                utility/libbacktrace/mmap.c \
                utility/libbacktrace/mmapio.c \
                utility/libbacktrace/posix.c \
                utility/libbacktrace/print.c \
                utility/libbacktrace/simple.c \
                utility/libbacktrace/sort.c \
                utility/libbacktrace/state.c \
                utility/libbacktrace/unknown.c
endif



SRC_INCDIRS  := . \
                utility

SUBMAKEFILES := 

ifeq ($(BUILDTESTS), 1)
SUBMAKEFILES += tests/bitsTest.mk \
                tests/filesTest.mk \
                tests/intervalListTest.mk \
                tests/loggingTest.mk \
                tests/stddevTest.mk
endif
