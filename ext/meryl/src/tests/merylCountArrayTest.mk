TARGET   := merylCountArrayTest
SOURCES  := merylCountArrayTest.C ../meryl/merylCountArray.C

SRC_INCDIRS  := . ../utility/src/utility ../meryl

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
