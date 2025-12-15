TARGET   := testHTSbam
SOURCES  := testHTSbam.C

SRC_INCDIRS := .. ../utility ../utility/htslib

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE} -lz -llzma -lbz2 -lcurl
TGT_PREREQS := lib${MODULE}.a
