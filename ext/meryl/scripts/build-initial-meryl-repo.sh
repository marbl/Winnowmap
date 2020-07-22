#!/bin/sh

#  You do not want to run this.  It was used to create the initial repo by
#  extracting the relevant bits from Canu/Meryl and renaming things.

if [ `pwd` != "/scratch/git/meryl-filtered" ] ; then
  echo Wrong directory.
  exit
fi

echo DELETE
rm -rf .git *

echo SYNC
rsync -a ../meryl-orig/ .

echo REWRITE
#../git-filter-repo/git-filter-repo --replace-refs delete-no-add \
#  --path-glob src/*.* \
#  --path src/meryl \
#  --path src/Makefile \
#  --path src/main.mk \
#  --path 'README.licenses' \
#  --path 'README.md'

#../git-filter-repo/git-filter-repo \
#  --path src/AS_UTL/Makefile \
#  --path src/AS_UTL/main.mk --invert-paths

../git-filter-repo/git-filter-repo \
  --path-glob 'src/AS_global.?' --invert-paths

../git-filter-repo/git-filter-repo \
  --path src/sequence --invert-paths

../git-filter-repo/git-filter-repo \
  --path-rename README.licenses:README.license

../git-filter-repo/git-filter-repo \
  --path-rename src/meryl_version_update.pl:scripts/version_update.pl

../git-filter-repo/git-filter-repo \
  --path-rename src/AS_MER:src/meryl

../git-filter-repo/git-filter-repo --path-rename       src/test-build.pl:src/tests/test-build.pl
../git-filter-repo/git-filter-repo --path-rename  src/test-operations.pl:src/tests/test-operations.pl

../git-filter-repo/git-filter-repo --path-rename src/meryl/meryl-check.C:src/meryl-check/meryl-check.C
../git-filter-repo/git-filter-repo --path-rename src/meryl/meryl-check.mk:src/meryl-check/meryl-check.mk

../git-filter-repo/git-filter-repo --path-rename src/meryl/meryl-import.C:src/meryl-import/meryl-import.C
../git-filter-repo/git-filter-repo --path-rename src/meryl/meryl-import.mk:src/meryl-import/meryl-import.mk

../git-filter-repo/git-filter-repo --path-rename src/meryl/meryl-lookup.C:src/meryl-lookup/meryl-lookup.C
../git-filter-repo/git-filter-repo --path-rename src/meryl/meryl-lookup.mk:src/meryl-lookup/meryl-lookup.mk


#../git-filter-repo/git-filter-repo --path-rename                     src/AS_global.H:src/utility/runtime.H
#../git-filter-repo/git-filter-repo --path-rename                     src/AS_global.C:src/utility/runtime.C
#../git-filter-repo/git-filter-repo --path-rename                     src/AS_global.h:src/utility/runtime.h
#../git-filter-repo/git-filter-repo --path-rename                     src/AS_global.c:src/utility/runtime.c

#  Delete all replace refs, since I can't seem to get them disabled.
git replace -d `git replace -l` > /dev/null 2>&1

exit

git remote add mmm /work/meryl
git fetch mmm --tags
git merge --allow-unrelated-histories mmm/master
git remote remove mmm

exit

Merges failed.


CONFLICT (add/add): Merge conflict in src/meryl/merylCountArray.C
Auto-merging src/meryl/merylCountArray.C
CONFLICT (add/add): Merge conflict in src/meryl/meryl-lookup.C
Auto-merging src/meryl/meryl-lookup.C
CONFLICT (add/add): Merge conflict in src/meryl/meryl-import.C
Auto-merging src/meryl/meryl-import.C
CONFLICT (add/add): Merge conflict in src/main.mk
Auto-merging src/main.mk
CONFLICT (add/add): Merge conflict in src/Makefile
Auto-merging src/Makefile
CONFLICT (add/add): Merge conflict in README.md
Auto-merging README.md
Recorded preimage for 'README.md'
Recorded preimage for 'src/Makefile'
Recorded preimage for 'src/main.mk'
Recorded preimage for 'src/meryl/meryl-import.C'
Recorded preimage for 'src/meryl/meryl-lookup.C'
Recorded preimage for 'src/meryl/merylCountArray.C'
Automatic merge failed; fix conflicts and then commit the result.

