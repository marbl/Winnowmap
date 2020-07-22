#!/bin/sh

#  You do not want to run this.  It was used to create the initial repo by
#  extracting the relevant bits from Canu and renaming things.

if [ `pwd` != "/scratch/git/meryl-rebuild" ] ; then
  echo Wrong directory.
  exit
fi

echo DELETE
rm -rf .git *

echo INIT
git init

echo MERGE CANU
git remote add ccc /scratch/git/canu-filtered
git fetch ccc --tags
git merge --allow-unrelated-histories ccc/master
git remote remove ccc

echo MERGE MERYL
git remote add mmm /scratch/git/meryl-filtered
git fetch mmm --tags
git merge --allow-unrelated-histories mmm/master
git remote remove mmm

cp -fp ../NEW-main.mk src/main.mk
cp -fp ../NEW-Makefile src/Makefile

cp -fp /scratch/git/canu-filtered/src/meryl/merylCountArray.C src/meryl/merylCountArray.C
cp -fp /scratch/git/canu-filtered/src/meryl/meryl-lookup.C src/meryl-lookup/meryl-lookup.C
cp -fp /scratch/git/canu-filtered/src/meryl/meryl-import.C src/meryl-import/meryl-import.C

git add -u
git status
git commit -m "Fix minor conflicts between meryl-repo and canu-repo."


git rm -r src/utility
git commit -m "Remove local copy of utility functions."

git submodule add git@github.com:marbl/meryl-utility src/utility
