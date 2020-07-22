#!/usr/bin/env perl

use strict;

my $nmeryl = "/work/meryl/FreeBSD-amd64/bin/meryl";
my $omeryl = "/work/canu/FreeBSD-amd64/bin/meryl";


sub makeSequence ($$$) {
    my $name   = shift @_;
    my $seq    = shift @_;
    my $seqLen = length($seq);
    my $lines  = shift @_;
    my $bases  = 0;

    if (-e $name) {
        for (my $ii=1; $ii <= $lines; $ii++) {
            $bases += $seqLen * $ii;
        }
    }

    else {
        open(F, "> $name");
        print F ">sequence\n";
        for (my $ii=1; $ii <= $lines; $ii++) {
            print F $seq x $ii;
            print F "\n";

            $bases += $seqLen * $ii;
        }
        close(F);
    }

    print STDERR "Generated '$name' with $bases bp of '$seq'.\n";

    return($bases);
}


if (0) {
    my $bases = makeSequence("A.fasta", "A", 5000);
    system("meryl k=22 print count n=$bases A.fasta output A > A.count.dump");
    system("meryl print A > A.dump");
}

if (0) {
    my $bases = makeSequence("AC.fasta", "AC", 5000);
    system("meryl k=22 print count n=$bases AC.fasta output AC > AC.count.dump");
    system("meryl print AC > AC.dump");
}

if (0) {
    my $bases = makeSequence("ACG.fasta", "ACG", 5000);
    system("meryl k=22 print count n=$bases ACG.fasta output ACG > ACG.count.dump");
    system("meryl print ACG > ACG.dump");
}

if (0) {
    my $bases = makeSequence("ACGT.fasta", "ACGT", 5000);
    system("meryl k=22 print count n=$bases ACGT.fasta output ACGT > ACGT.count.dump");
    system("meryl print ACGT > ACGT.dump");
}



system("$omeryl -B -C -m 22 -s DATA/reads.fasta -o old-reads");
system("$omeryl -Dt -s old-reads > old-reads.dump");

foreach my $memory (qw(01 02 04 08 16 32)) {
foreach my $thread (qw(01 12)) {
    system("$nmeryl k=22 threads=$thread memory=$memory count DATA/reads.fasta output new-reads-mem$memory-thr$thread");
    system("$nmeryl print new-reads-mem$memory-thr$thread > new-reads-mem$memory-thr$thread.dump");
}
}
