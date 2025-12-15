#!/bin/sh
#
#  Before building a release:
#
#  Make a place to work, grab the bits you want to release:
#    git clone git@github.com:marbl/PACKAGE PACKAGE-release
#    cd PACKAGE-release
#
#  Commit to master:
#    Increase version in documentation/source/conf.py  (if present)
#    Increase version in src/main.mk
#
#  Pull in submodule code (as a side-effect of 'clean'):
#    cd src && gmake clean
#
#  Tag the next release development at the tip in master:
#    git tag -a v2.2-development -m "Development for v2.2."
#    git push --follow-tags
#
#  Make a branch:
#    git checkout -b v2.1-maintenance
#
#  Commit to branch:
#    Change 'snapshot' to 'release' in src/main.mk (VERSION := release v2.1)
#    git push --set-upstream origin v2.1-maintenance
#
#  Run this script:
#    scripts/buildRelease.sh PACKAGE 2.1
#
#  Some notes:
#   - The source code from the branch made above will be used
#     to build tarballs with source and binaries.
#   - These tarballs (and scratch space) are called 'build-*'.
#   - The very first step is to rename '.git' to 'dot-git-directory';
#     catastrophic failure (of this script) could leave it
#     renamed.
#
#  Ubuntu needs zlib1g-dev, libcurl4-openssl-dev, libssl-dev, liblzma-dev,
#  libbz2-dev.
#

set -e
#et -x

if [ -e src/main.mk ] ; then
  package=$( grep -e '^MODULE[[:space:]]*:='  src/main.mk | awk '{ print $3 }' )
  release=$( grep -e '^VERSION[[:space:]]*:=' src/main.mk | awk '{ print $3 }' )
  version=$( grep -e '^VERSION[[:space:]]*:=' src/main.mk | awk '{ print $4 }' | sed s/v// )
else                                  echo "ERROR: no src/main.mk found."                                         ; exit 1 ; fi

if [ "$package"  = "" ]        ; then echo "ERROR: no MODULE name found in src/main.mk."                          ; exit 1 ; fi
if [ "$version"  = "" ]        ; then echo "ERROR: no VERSION number found in src/main.mk."                       ; exit 1 ; fi
if [ "$release" != "release" ] ; then echo "ERROR: expecting 'release' VERSION in src/main.mk, found '$release'." ; exit 1 ; fi

#
#
echo "Building release for '$package' version '$version'."

if [ -e .git ] ; then
  echo "Moving .git directory out of the way."
  mv .git dot-git-directory
fi

#
#
echo "Preparing build trees."

rm -rf build
rm -rf build-darwin build-darwin.out
rm -rf build-linux  build-linux.out
rm -rf build-src

rm  -f build-linux.sh

rm  -f ${package}-${version}.Darwin-aarch64.tar ${package}-${version}.Darwin-aarch64.tar.xz
rm  -f ${package}-${version}.Linux-amd64.tar    ${package}-${version}.Linux-amd64.tar.xz
rm  -f ${package}-${version}.tar                ${package}-${version}.tar.xz

mkdir -p build-src    ; rsync -a src/ build-src/src    ; cp -p README* build-src/
mkdir -p build-darwin ; rsync -a src/ build-darwin/src ; cp -p README* build-darwin/
mkdir -p build-linux  ; rsync -a src/ build-linux/src  ; cp -p README* build-linux/


#
#  HBB needs bzip2 (only for libbz2 though) and xz (for liblzma) installed
#  manually.  Download those.
#
#  This takes anywhere from 2 1/2 to 10 minutes (!) on my M1 Air (in Docker,
#  running as amd64) depending on if there is _ANYTHING_ else running on the
#  machine.
#
if [ ! -e bzip2-1.0.8.tar.gz ] ; then
  curl -LRO https://sourceware.org/pub/bzip2/bzip2-1.0.8.tar.gz
fi
if [ ! -e bzip2-1.0.8 ] ; then
  gzip -dc bzip2-1.0.8.tar.gz | tar -xf -
fi

if [ ! -e xz-5.6.3.tar.gz ] ; then
  curl -LRO https://github.com/tukaani-project/xz/releases/download/v5.6.3/xz-5.6.3.tar.gz
fi
if [ ! -e xz-5.6.3 ] ; then
  gzip -dc xz-5.6.3.tar.gz | tar -xf -
fi


#
#
echo >> build-linux.sh  "#!/bin/bash"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "set -e"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "#  HBB docs say that this will inflate binary sizes, and instead we"
echo >> build-linux.sh  "#  need to install from source."
echo >> build-linux.sh  "#  Didn't seem to matter, given that debug symbols are GIGANTIC now."
echo >> build-linux.sh  "#yum install -y xz-devel bzip2-devel"
echo >> build-linux.sh  ""
echo >> build-linux.sh  ""
echo >> build-linux.sh  "if [ -e /dock/hbb-xz-bz2-install.tar ] ; then"
echo >> build-linux.sh  "  echo Installing libbz2.a and liblzma.a."
echo >> build-linux.sh  ""
echo >> build-linux.sh  "  cd /hbb_exe"
echo >> build-linux.sh  "  tar -xf /dock/hbb-xz-bz2-install.tar"
echo >> build-linux.sh  "  cd /dock"
echo >> build-linux.sh  "fi"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "#  This is pretty quick..."
echo >> build-linux.sh  "if [ ! -e /hbb_exe/lib/libbz2.a ] ; then"
echo >> build-linux.sh  "  echo Building libbz2.a."
echo >> build-linux.sh  ""
echo >> build-linux.sh  "  cd /dock/bzip2-1.0.8"
echo >> build-linux.sh  "  gcc $STATICLIB_CFLAGS -c -o blocksort.o  blocksort.c"
echo >> build-linux.sh  "  gcc $STATICLIB_CFLAGS -c -o huffman.o    huffman.c"
echo >> build-linux.sh  "  gcc $STATICLIB_CFLAGS -c -o crctable.o   crctable.c"
echo >> build-linux.sh  "  gcc $STATICLIB_CFLAGS -c -o randtable.o  randtable.c"
echo >> build-linux.sh  "  gcc $STATICLIB_CFLAGS -c -o compress.o   compress.c"
echo >> build-linux.sh  "  gcc $STATICLIB_CFLAGS -c -o decompress.o decompress.c"
echo >> build-linux.sh  "  gcc $STATICLIB_CFLAGS -c -o bzlib.o      bzlib.c"
echo >> build-linux.sh  "  ar cq libbz2.a blocksort.o huffman.o crctable.o randtable.o compress.o decompress.o bzlib.o"
echo >> build-linux.sh  "  cp -f libbz2.a /hbb_exe/lib/"
echo >> build-linux.sh  "  cp -f bzlib.h  /hbb_exe/include"
echo >> build-linux.sh  "fi"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "#  ...but this can take up to 8 minutes if the machine isn't idle."
echo >> build-linux.sh  "if [ ! -e /hbb_exe/lib/liblzma.a ] ; then"
echo >> build-linux.sh  "  echo Building liblzma.a.  This is very slow."
echo >> build-linux.sh  ""
echo >> build-linux.sh  "  export CFLAGSPREV=$CFLAGS"
echo >> build-linux.sh  "  export CFLAGS=$STATICLIB_CFLAGS"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "  cd /dock/xz-5.6.3"
echo >> build-linux.sh  "  ./configure \\"
echo >> build-linux.sh  "    --prefix=/hbb_exe \\"
echo >> build-linux.sh  "    --disable-shared --enable-static \\"
echo >> build-linux.sh  "    --disable-xz --disable-xzdec --disable-lzmadec --disable-lzmainfo \\"
echo >> build-linux.sh  "    --disable-lzma-links --disable-scripts --disable-doc --disable-doxygen --disable-nls > configure.out 2>&1"
echo >> build-linux.sh  "  gmake --no-print-directory -j 4 install > gmake.out 2>&1"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "  export CFLAGS=$CFLAGSPREV"
echo >> build-linux.sh  "  unset  CFLAGSPREV"
echo >> build-linux.sh  "fi"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "#  Save the work so we can reuse it on reruns."
echo >> build-linux.sh  "if [ ! -e /dock/hbb-xz-bz2-install.tar ] ; then"
echo >> build-linux.sh  "  echo Saving libbz2.a and liblzma.a installations for later reuse."
echo >> build-linux.sh  ""
echo >> build-linux.sh  "  cd /hbb_exe"
echo >> build-linux.sh  "  tar -cf /dock/hbb-xz-bz2-install.tar \\"
echo >> build-linux.sh  "    lib/libbz2.a \\"
echo >> build-linux.sh  "    include/bzlib.h \\"
echo >> build-linux.sh  "    lib/liblzma* \\"
echo >> build-linux.sh  "    lib/pkgconfig/liblzma.pc \\"
echo >> build-linux.sh  "    include/lzma*"
echo >> build-linux.sh  "  cd /dock"
echo >> build-linux.sh  "fi"
echo >> build-linux.sh  ""
echo >> build-linux.sh  ""
echo >> build-linux.sh  "#  Fix for HBB's including debug symbols in libraries we don't care about."
echo >> build-linux.sh  "cd /hbb_exe/lib"
echo >> build-linux.sh  "strip -d libbz2.a"
echo >> build-linux.sh  "strip -d libcurl.a"
echo >> build-linux.sh  "strip -d liblzma.a"
echo >> build-linux.sh  "strip -d libsqlite3.a"
echo >> build-linux.sh  "strip -d libstdc++.a"
echo >> build-linux.sh  "strip -d libsupc++.a"
echo >> build-linux.sh  "strip -d libz.a"
echo >> build-linux.sh  "ranlib *.a"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "cd /hbb_exe/lib64"
echo >> build-linux.sh  "strip -d libcrypto.a"
echo >> build-linux.sh  "strip -d libssl.a"
echo >> build-linux.sh  "ranlib *.a"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "#  Fix for HBB's gcc 9.2.1 mishandling of scoped enums and bitfields."
echo >> build-linux.sh  "#  https://stackoverflow.com/questions/36005063/gcc-suppress-warning-too-small-to-hold-all-values-of"
echo >> build-linux.sh  "export CXXFLAGS=\"\$CXXFLAGS -fpermissive\""
echo >> build-linux.sh  ""
echo >> build-linux.sh  ""
echo >> build-linux.sh  "rm -rf /dock/build"
echo >> build-linux.sh  ""
echo >> build-linux.sh  ""
echo >> build-linux.sh  "#  Build htslib and parasail WITHOUT debug symbols."
echo >> build-linux.sh  "#  There doesn't appear to be an easy way to do this, and"
echo >> build-linux.sh  "#  we just build libcanu.a fully, then remove the stuff"
echo >> build-linux.sh  "#  we want to recompile with symbols."
echo >> build-linux.sh  "#"
echo >> build-linux.sh  "echo Building ${package} \(without debug\)."
echo >> build-linux.sh  "cd /dock/src"
echo >> build-linux.sh  "gmake -j 4 BUILDOPTIMIZED=1 /dock/build/lib/libcanu.a > ../build-linux-debug.out 2>&1"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "cd /dock/build/obj/lib/lib*.a"
echo >> build-linux.sh  "rm -rf \`ls -d *             | grep -Ev utility\`"
echo >> build-linux.sh  "rm -rf \`ls -d utility/src/* | grep -Ev htslib\|parasail\`"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "rm -f /dock/build/lib/lib*.a"
echo >> build-linux.sh  ""
echo >> build-linux.sh  ""
echo >> build-linux.sh  "#  Recompile the rest with debug sybols."
echo >> build-linux.sh  "#"
echo >> build-linux.sh  "echo Building ${package}."
echo >> build-linux.sh  "cd /dock/src"
echo >> build-linux.sh  "gmake -j 4 > ../build-linux.out 2>&1"
echo >> build-linux.sh  "cd .."
echo >> build-linux.sh  ""
echo >> build-linux.sh  "mv build/* build-linux/"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "rm -rf build-darwin/obj"
echo >> build-linux.sh  "rm -rf build-linux/obj"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "mv build-darwin ${package}-${version}"
echo >> build-linux.sh  "tar -cf ${package}-${version}.Darwin-aarch64.tar ${package}-${version}/README* ${package}-${version}/bin ${package}-${version}/lib ${package}-${version}/share"
echo >> build-linux.sh  "mv ${package}-${version} build-darwin"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "mv build-linux ${package}-${version}"
echo >> build-linux.sh  "tar -cf ${package}-${version}.Linux-amd64.tar  ${package}-${version}/README*  ${package}-${version}/bin  ${package}-${version}/lib  ${package}-${version}/share"
echo >> build-linux.sh  "mv ${package}-${version} build-linux"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "mv build-src ${package}-${version}"
echo >> build-linux.sh  "tar -cf ${package}-${version}.tar              ${package}-${version}/README*  ${package}-${version}/src"
echo >> build-linux.sh  "mv ${package}-${version} build-src"
echo >> build-linux.sh  ""
echo >> build-linux.sh  ""

chmod 755 build-linux.sh

#
#
echo "Build for MacOS."

cd src
gmake -j 4 BUILDRELEASE=1 > ../build-darwin.out 2>&1
cd ..

mv build/* build-darwin/

#
#
echo "Build for Linux and make tarballs."

echo \
docker run --platform=linux/amd64 -v `pwd`:/dock -t -i --rm phusion/holy-build-box:latest /hbb_exe/activate-exec bash /dock/build-linux.sh
docker run --platform=linux/amd64 -v `pwd`:/dock -t -i --rm phusion/holy-build-box:latest /hbb_exe/activate-exec bash /dock/build-linux.sh

#
#
echo "Compress."

xz -9v ${package}-${version}.Darwin-aarch64.tar
xz -9v ${package}-${version}.Linux-amd64.tar
xz -9v ${package}-${version}.tar

if [ -e dot-git-directory ] ; then
  echo "Restoring .git directory."
  mv dot-git-directory .git
fi

#
#
echo "Finished."
exit
