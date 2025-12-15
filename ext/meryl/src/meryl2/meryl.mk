TARGET   := meryl2
SOURCES  := meryl.C \
            merylCommandBuilder-isAssign.C \
            merylCommandBuilder-isSelect.C \
            merylCommandBuilder-printTree.C \
            merylCommandBuilder-processText.C \
            merylCommandBuilder-spawnThreads.C \
            merylCommandBuilder.C \
            merylCountArray.C \
            merylGlobals.C \
            merylInput.C \
            merylInput-canu.C \
            merylSelector.C \
            merylOp-count-memorySize.C \
            merylOp-count.C \
            merylOp-countSequential.C \
            merylOp-countSimple.C \
            merylOp-countThreads.C \
            merylOp-nextMer.C \
            merylOp.C \
            merylOpCompute.C \
            merylOpTemplate.C

SRC_INCDIRS := .

#  If we're part of Canu, build with canu support and use Canu's copy of
#  meryl-utility.  Otherwise, don't.
ifneq ($(wildcard stores/sqStore.H), )
  SRC_CXXFLAGS := -DCANU
  SRC_INCDIRS  := ../../../utility/src ../../../stores

#  If we're part of something else, include the something else's
#  utility directory.
else ifneq ($(wildcard meryl/src/meryl2/meryl.C), )
  SRC_INCDIRS  := ../../../utility/src

#  Otherwise, we're building directly in the meryl repo.
else
  SRC_INCDIRS  := ../utility/src

endif


TGT_LDFLAGS  := -L${TARGET_DIR}/lib
TGT_LDLIBS   := -l${MODULE}
TGT_PREREQS  := lib${MODULE}.a
