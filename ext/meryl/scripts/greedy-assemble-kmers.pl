#!/usr/bin/env perl

################################################################################
 #
 #  This file is part of meryl, a genomic k-kmer counter with nice features.
 #
 #  This software is based on:
 #    'Canu' v2.0              (https://github.com/marbl/canu)
 #  which is based on:
 #    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 #
 #  Except as indicated otherwise, this is a 'United States Government Work',
 #  and is released in the public domain.
 #
 #  File 'README.licenses' in the root directory of this distribution
 #  contains full conditions and disclaimers.
 ##

use strict;

my %kmers;
my $slen;

#  Reverse complement a string.
#
sub rc ($) {
    my $s = $_[0];

    $s =~ tr/ACTG/TGAC/;
    $s = reverse $s;

    return($s);
}


#  Load kmers into memory.
#
sub loadKmers () {
    while (<STDIN>) {
        my ($k, $c) = split '\s+', $_;

        $kmers{$k} = 1;
    }

    $slen = length((keys %kmers)[0]) - 1;
}


#  Pick a random kmer to be the seed of a contig, and remove it from the set
#  of kmers to place.
#
sub pickSeed () {
    my $seed = (keys %kmers)[0];

    delete $kmers{$seed};

    return($seed);
}


#  Pull k-1 bases from each end.  We'll try to extend each by a single
#  base to greedily grow the contig.
#
sub extend ($) {
    my $asm = shift @_;

    my $l = substr($asm,      0, $slen);
    my $r = substr($asm, -$slen, $slen);

    foreach my $e ("A", "C", "T", "G") {
        my $lkf = "$e$l";
        my $lkr = rc($lkf);

        if (exists($kmers{$lkf}))   {  delete $kmers{$lkf};  return(1, "$e$asm");  }
        if (exists($kmers{$lkr}))   {  delete $kmers{$lkr};  return(1, "$e$asm");  }
    }

    foreach my $e ("A", "C", "T", "G") {
        my $rkf = "$r$e";
        my $rkr = rc($rkf);

        if (exists($kmers{$rkf}))   {  delete $kmers{$rkf};  return(1, "$asm$e");  }
        if (exists($kmers{$rkr}))   {  delete $kmers{$rkr};  return(1, "$asm$e");  }
    }

    return(0, $asm);
}


#
#  MAIN!
#

loadKmers();

my $idx = 1;

while (scalar(keys %kmers) > 0) {
    my $asm = pickSeed();

    while (1) {
        my $sf;
        my $sr;

        ($sf, $asm) = extend($asm);
        $asm = rc($asm);

        ($sr, $asm) = extend($asm);
        $asm = rc($asm);

        last   if ($sf + $sr == 0);
    }

    print ">$idx\n";  $idx++;
    print "$asm\n";
}

exit(0);
