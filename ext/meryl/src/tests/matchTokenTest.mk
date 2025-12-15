TARGET   := matchTokenTest
SOURCES  := matchTokenTest.C

SRC_INCDIRS  := . ../utility/src ../meryl2

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
