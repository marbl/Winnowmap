TARGET   := testVectorSupport
SOURCES  := testVectorSupport.C

SRC_INCDIRS := .. ../utility

TGT_CXXFLAGS:= -mxsave
TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
