#!/bin/sh

mkdir -p test
mkdir -p test/faidx

cat > test/emptyfile <<EOF
EOF

cat > test/xx#pair.sam <<EOF
@SQ     SN:xx   LN:20
a1      99      xx      1       1       10M     =       11      20      AAAAAAAAAA      **********
b1      99      xx      1       1       10M     =       11      20      AAAAAAAAAA      **********
c1      99      xx      1       1       10M     =       11      20      AAAAAAAAAA      **********
a1      147     xx      11      1       10M     =       1       -20     TTTTTTTTTT      **********
b1      147     xx      11      1       10M     =       1       -20     TTTTTTTTTT      **********
c1      147     xx      11      1       10M     =       1       -20     TTTTTTTTTT      **********
EOF


for f in "emptyfile" "ce.fa" "xx.fa" "faidx/fastqs.fq" ; do
  if [ ! -e test/$f ] ; then
    curl -LR "https://raw.githubusercontent.com/samtools/htslib/develop/test/$f" -o "test/$f"
  fi
done

../../build/bin/testHTSlib && echo Pass\!
