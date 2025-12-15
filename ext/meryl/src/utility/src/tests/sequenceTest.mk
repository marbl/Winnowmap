TARGET   := sequenceTest
SOURCES  := sequenceTest.C

SRC_INCDIRS := .. ../utility

#  htslib needs the compressors and curl.
TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE} -lz -llzma -lbz2 -lcurl
TGT_PREREQS := lib${MODULE}.a
