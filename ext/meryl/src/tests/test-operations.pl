#!/usr/bin/env perl

use strict;

my $nmeryl = "/work/meryl/FreeBSD-amd64/bin/meryl";
my $omeryl = "/work/canu/FreeBSD-amd64/bin/meryl";

testBuild();
testUnion();
testIntersect();
exit(0);



sub testBuild {

    system("mkdir new");

    system("$nmeryl k=22 count n=5000000 DATA/ecoli-dh.fasta output new/ecoli-dh");
    system("$nmeryl k=22 count n=5000000 DATA/ecoli-mg.fasta output new/ecoli-mg");
    system("$nmeryl k=22 count n=5000000 DATA/ecoli-w3.fasta output new/ecoli-w3");

    system("$nmeryl print new/ecoli-dh > new/ecoli-dh.dump");
    system("$nmeryl print new/ecoli-mg > new/ecoli-mg.dump");
    system("$nmeryl print new/ecoli-w3 > new/ecoli-w3.dump");

    system("mkdir old");

    system("$omeryl -B -C -m 22 -s DATA/ecoli-dh.fasta -o old/ecoli-dh");
    system("$omeryl -B -C -m 22 -s DATA/ecoli-mg.fasta -o old/ecoli-mg");
    system("$omeryl -B -C -m 22 -s DATA/ecoli-w3.fasta -o old/ecoli-w3");

    system("$omeryl -Dt -s old/ecoli-dh >  old/ecoli-dh.dump");
    system("$omeryl -Dt -s old/ecoli-mg >  old/ecoli-mg.dump");
    system("$omeryl -Dt -s old/ecoli-w3 >  old/ecoli-w3.dump");

    my %ndump;

    open(F, "md5 -r new/*dump |");
    while(<F>) {
        if (m/^(\w+)\s+new\/(.*.dump)$/) {
            $ndump{$2} = $1;
        }
    }
    close(F);

    my $fail = 0;

    open(F, "md5 -r old/*dump |");
    while(<F>) {
        if (m/^(\w+)\s+old\/(.*.dump)$/) {
            if ($ndump{$2} ne $1) {
                print "$2 dump md5 fails: new $ndump{$2} != old $1\n";
                $fail++;
            }
        }
    }
    close(F);

    die "testBuild() fails.\n"  if ($fail > 0);

    print STDERR "testBuild() PASSES!\n";
}



sub testUnion {
    my ($a, $am, $ac);
    my ($b, $bm, $bc);
    my ($c, $cm, $cc);

    print STDERR "Generating UNION with meryl.\n";

    system("meryl union     new/ecoli-dh new/ecoli-mg new/ecoli-w3 output new/union");
    system("meryl union-min new/ecoli-dh new/ecoli-mg new/ecoli-w3 output new/union-min");
    system("meryl union-max new/ecoli-dh new/ecoli-mg new/ecoli-w3 output new/union-max");
    system("meryl union-sum new/ecoli-dh new/ecoli-mg new/ecoli-w3 output new/union-sum");

    open(NEXI, "meryl print new/union     |");
    open(NMIN, "meryl print new/union-min |");
    open(NMAX, "meryl print new/union-max |");
    open(NSUM, "meryl print new/union-sum |");

    print STDERR "Comparing UNION against stuipid algorithm.\n";

    open(A, "< new/ecoli-dh.dump") or die;
    open(B, "< new/ecoli-mg.dump") or die;
    open(C, "< new/ecoli-w3.dump") or die;

    #open(UEXI, "> eval-union")     or die;
    #open(UMIN, "> eval-union-min") or die;
    #open(UMAX, "> eval-union-max") or die;
    #open(USUM, "> eval-union-sum") or die;

    while (!eof(A) || !eof(B) || !eof(C)) {

        #  Load the next mer, if needed.

        if (!defined($a)) {
            $a = <A>;
            ($am, $ac) = split '\s+', $a;
        }

        if (!defined($b)) {
            $b = <B>;
            ($bm, $bc) = split '\s+', $b;
        }

        if (!defined($c)) {
            $c = <C>;
            ($cm, $cc) = split '\s+', $c;
        }

        #  To make perl's lexicographic comparison correct, we need to swap G and T.

        $am =~ tr/GT/TG/;
        $bm =~ tr/GT/TG/;
        $cm =~ tr/GT/TG/;

        #  Find the smallest kmer.

        my $mm;

        $mm = $am  if (defined($am));
        $mm = $bm  if (defined($bm) && ($bm le $mm));
        $mm = $cm  if (defined($cm) && ($cm le $mm));

        #
        #  UNION!
        #

        my $cexi = 0;
        my $cmin = 999999999;
        my $cmax = 0;
        my $csum = 0;

        if ($am eq $mm) {
            $cexi  = 1;
            $cmin  = $ac  if ($ac < $cmin);
            $cmax  = $ac  if ($cmax < $ac);
            $csum += $ac;
        }

        if ($bm eq $mm) {
            $cexi  = 1;
            $cmin  = $bc  if ($bc < $cmin);
            $cmax  = $bc  if ($cmax < $bc);
            $csum += $bc;
        }

        if ($cm eq $mm) {
            $cexi  = 1;
            $cmin  = $cc  if ($cc < $cmin);
            $cmax  = $cc  if ($cmax < $cc);
            $csum += $cc;
        }


        #  Undo the letter swap.

        $am =~ tr/GT/TG/;
        $bm =~ tr/GT/TG/;
        $cm =~ tr/GT/TG/;

        $mm =~ tr/GT/TG/;

        #  Dump our version.

        #print UEXI "$mm\t$cexi\n";
        #print UMIN "$mm\t$cmin\n";
        #print UMAX "$mm\t$cmax\n";
        #print USUM "$mm\t$csum\n";

        #  Compare against what meryl said.

        if ($cexi > 0) {
            my $nexi = <NEXI>;  chomp $nexi;
            my $nmin = <NMIN>;  chomp $nmin;
            my $nmax = <NMAX>;  chomp $nmax;
            my $nsum = <NSUM>;  chomp $nsum;

            my $fail = 0;

            $fail++   if ($nexi ne "$mm\t$cexi");
            $fail++   if ($nmin ne "$mm\t$cmin");
            $fail++   if ($nmax ne "$mm\t$cmax");
            $fail++   if ($nsum ne "$mm\t$csum");

            if ($fail > 0) {
                print STDERR "FAILED on mm=$mm\n";
                print STDERR "FAILED\n";
                print STDERR "FAILED on am=$am $ac\n";
                print STDERR "FAILED on bm=$bm $bc\n";
                print STDERR "FAILED on cm=$cm $cc\n";
                print STDERR "FAILED\n";
                print STDERR "FAILED nexi '$nexi'\n";
                print STDERR "FAILED nmin '$nmin'\n";
                print STDERR "FAILED nmax '$nmax'\n";
                print STDERR "FAILED nsum '$nsum'\n";
                exit(1);
            }
        }

        #  Setup for the next kmer

        if ($am eq $mm) { $a = undef; $am = undef; $ac = 0; }
        if ($bm eq $mm) { $b = undef; $bm = undef; $bc = 0; }
        if ($cm eq $mm) { $c = undef; $cm = undef; $cc = 0; }
    }

    close(USUM);
    close(UMAX);
    close(UMIN);
    close(UEXI);

    close(C);
    close(B);
    close(A);

    print STDERR "testUnion() PASSES!\n";
}



sub testIntersect {
    my ($a, $am, $ac);
    my ($b, $bm, $bc);
    my ($c, $cm, $cc);

    print STDERR "Generating INTERSECT with meryl.\n";

    system("meryl intersect     new/ecoli-dh new/ecoli-mg new/ecoli-w3 output new/intersect");
    system("meryl intersect-min new/ecoli-dh new/ecoli-mg new/ecoli-w3 output new/intersect-min");
    system("meryl intersect-max new/ecoli-dh new/ecoli-mg new/ecoli-w3 output new/intersect-max");
    system("meryl intersect-sum new/ecoli-dh new/ecoli-mg new/ecoli-w3 output new/intersect-sum");

    open(NEXI, "meryl print new/intersect     |");
    open(NMIN, "meryl print new/intersect-min |");
    open(NMAX, "meryl print new/intersect-max |");
    open(NSUM, "meryl print new/intersect-sum |");

    print STDERR "Comparing INTERSECT against stuipid algorithm.\n";

    open(A, "< new/ecoli-dh.dump") or die;
    open(B, "< new/ecoli-mg.dump") or die;
    open(C, "< new/ecoli-w3.dump") or die;

    #open(UEXI, "> eval-intersect")     or die;
    #open(UMIN, "> eval-intersect-min") or die;
    #open(UMAX, "> eval-intersect-max") or die;
    #open(USUM, "> eval-intersect-sum") or die;

    while (!eof(A) || !eof(B) || !eof(C)) {

        #  Load the next mer, if needed.

        if (!defined($a)) {
            $a = <A>;
            ($am, $ac) = split '\s+', $a;
        }

        if (!defined($b)) {
            $b = <B>;
            ($bm, $bc) = split '\s+', $b;
        }

        if (!defined($c)) {
            $c = <C>;
            ($cm, $cc) = split '\s+', $c;
        }

        #  To make perl's lexicographic comparison correct, we need to swap G and T.

        $am =~ tr/GT/TG/;
        $bm =~ tr/GT/TG/;
        $cm =~ tr/GT/TG/;

        #  Find the smallest kmer.

        my $mm;

        $mm = $am  if (defined($am));
        $mm = $bm  if (defined($bm) && ($bm le $mm));
        $mm = $cm  if (defined($cm) && ($cm le $mm));

        #
        #  INTERSECT!
        #

        my $cexi = 0;
        my $cmin = 999999999;
        my $cmax = 0;
        my $csum = 0;

        if (($mm eq $am) &&
            ($mm eq $bm) &&
            ($mm eq $cm)) {
            $cexi  = 1;

            $cmin  = $ac  if ($ac < $cmin);
            $cmin  = $bc  if ($bc < $cmin);
            $cmin  = $cc  if ($cc < $cmin);

            $cmax  = $ac  if ($cmax < $ac);
            $cmax  = $bc  if ($cmax < $bc);
            $cmax  = $cc  if ($cmax < $cc);

            $csum  = $ac + $bc + $cc;
        }
 
        #  Undo the letter swap.

        $am =~ tr/GT/TG/;
        $bm =~ tr/GT/TG/;
        $cm =~ tr/GT/TG/;

        $mm =~ tr/GT/TG/;

        #  Dump our version.

        #print UEXI "$mm\t$cexi\n";
        #print UMIN "$mm\t$cmin\n";
        #print UMAX "$mm\t$cmax\n";
        #print USUM "$mm\t$csum\n";

        #  Compare against what meryl said.

        if ($cexi > 0) {
            my $nexi = <NEXI>;  chomp $nexi;
            my $nmin = <NMIN>;  chomp $nmin;
            my $nmax = <NMAX>;  chomp $nmax;
            my $nsum = <NSUM>;  chomp $nsum;

            my $fail = 0;

            $fail++   if ($nexi ne "$mm\t$cexi");
            $fail++   if ($nmin ne "$mm\t$cmin");
            $fail++   if ($nmax ne "$mm\t$cmax");
            $fail++   if ($nsum ne "$mm\t$csum");

            if ($fail > 0) {
                print STDERR "FAILED on mm=$mm\n";
                print STDERR "FAILED\n";
                print STDERR "FAILED on am=$am $ac\n";
                print STDERR "FAILED on bm=$bm $bc\n";
                print STDERR "FAILED on cm=$cm $cc\n";
                print STDERR "FAILED\n";
                print STDERR "FAILED nexi '$nexi'\n";
                print STDERR "FAILED nmin '$nmin'\n";
                print STDERR "FAILED nmax '$nmax'\n";
                print STDERR "FAILED nsum '$nsum'\n";
            }
        }

        #  Setup for the next kmer

        if ($am eq $mm) { $a = undef; $am = undef; $ac = 0; }
        if ($bm eq $mm) { $b = undef; $bm = undef; $bc = 0; }
        if ($cm eq $mm) { $c = undef; $cm = undef; $cc = 0; }
    }

    close(USUM);
    close(UMAX);
    close(UMIN);
    close(UEXI);

    close(C);
    close(B);
    close(A);

    print STDERR "testIntersect() PASSES!\n";
}
