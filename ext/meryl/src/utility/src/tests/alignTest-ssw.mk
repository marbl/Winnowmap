TARGET   := alignTest-ssw
SOURCES  := alignTest-ssw.C

SRC_INCDIRS := .. ../utility

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
