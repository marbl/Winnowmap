TARGET   := parasailTest
SOURCES  := parasailTest.C

#SRC_INCDIRS := .. ../utility /work/software/parasail
SRC_INCDIRS := .. ../utility ../parasail

#TGT_LDFLAGS := -L${TARGET_DIR}/lib -L/work/software/parasail/.libs -lparasail
TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
