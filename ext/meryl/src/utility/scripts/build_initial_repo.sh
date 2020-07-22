#!/bin/sh

#  You do not want to run this.  It was used to create the initial repo by
#  extracting the relevant bits from Canu and renaming things.

if [ `pwd` != "/scratch/git/canu-test" ] ; then
  echo Wrong directory.
  exit
fi

echo DELETE
rm -rf .git *

echo SYNC
rsync -a ../canu-orig/ .

echo REWRITE
../git-filter-repo/git-filter-repo --replace-refs delete-no-add \
  --path kmer/libutil \
  --path kmer/libkmer \
  --path kmer/libbio \
  --path src/AS_UTL \
  --path src/AS_global.C \
  --path src/AS_global.H \
  --path src/AS_global.c \
  --path src/AS_global.h \
  --path src/utility \
  --path src/canu_version_update.pl \
  --path src/Makefile \
  --path src/main.mk \
  --path 'README.license.GPL' \
  --path 'README.licenses' \
  --path 'README.md'

../git-filter-repo/git-filter-repo \
  --path src/AS_UTL/Makefile \
  --path src/AS_UTL/main.mk --invert-paths

../git-filter-repo/git-filter-repo \
  --path-rename src/canu_version_update.pl:scripts/version_update.pl

../git-filter-repo/git-filter-repo \
  --path-rename src/AS_UTL:src/utility

../git-filter-repo/git-filter-repo --path-rename                     src/AS_global.H:src/utility/runtime.H
../git-filter-repo/git-filter-repo --path-rename                     src/AS_global.C:src/utility/runtime.C
../git-filter-repo/git-filter-repo --path-rename                     src/AS_global.h:src/utility/runtime.h
../git-filter-repo/git-filter-repo --path-rename                     src/AS_global.c:src/utility/runtime.c

../git-filter-repo/git-filter-repo --path-rename              src/utility/bitsTest.C:src/tests/bitsTest.C
../git-filter-repo/git-filter-repo --path-rename             src/utility/bitsTest.mk:src/tests/bitsTest.mk
../git-filter-repo/git-filter-repo --path-rename             src/utility/filesTest.C:src/tests/filesTest.C
../git-filter-repo/git-filter-repo --path-rename            src/utility/filesTest.mk:src/tests/filesTest.mk
../git-filter-repo/git-filter-repo --path-rename      src/utility/intervalListTest.C:src/tests/intervalListTest.C
../git-filter-repo/git-filter-repo --path-rename     src/utility/intervalListTest.mk:src/tests/intervalListTest.mk
../git-filter-repo/git-filter-repo --path-rename           src/utility/loggingTest.C:src/tests/loggingTest.C
../git-filter-repo/git-filter-repo --path-rename          src/utility/loggingTest.mk:src/tests/loggingTest.mk
../git-filter-repo/git-filter-repo --path-rename  src/utility/memoryMappedFileTest.C:src/tests/memoryMappedFileTest.C
../git-filter-repo/git-filter-repo --path-rename         src/utility/mt19937arTest.C:src/tests/mt19937arTest.C
../git-filter-repo/git-filter-repo --path-rename            src/utility/stddevTest.C:src/tests/stddevTest.C
../git-filter-repo/git-filter-repo --path-rename           src/utility/stddevTest.mk:src/tests/stddevTest.mk
../git-filter-repo/git-filter-repo --path-rename           src/utility/stringsTest.C:src/tests/stringsTest.C
../git-filter-repo/git-filter-repo --path-rename              src/utility/testRand.C:src/tests/randTest.C


#  Delete all replace refs, since I can't seem to get them disabled.
git replace -d `git replace -l` > refs-deleted 2>&1


cp -fp ../NEW-Makefile          src/Makefile
cp -fp ../NEW-main.mk           src/main.mk
cp -fp ../NEW-runtime.C         src/utility/runtime.C
cp -fp ../NEW-runtime.H         src/utility/runtime.H
cp -fp ../NEW-types.H           src/utility/types.H
cp -fp ../NEW-version_update.pl scripts/version_update.pl
