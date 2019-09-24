#!/usr/bin/env perl
# TRF_blast_group.pl
# use grouped tandem fasta to identify regions in the reduced.trf.txt file

use warnings;
use strict;
use POSIX;
use List::Util qw( min max );
use List::Util qw(sum);
use Math::Round;

open FA, "$ARGV[0]" or die; ## .fa of the grouped TRF file.
open TRF, "$ARGV[1]" or die; ## processed TRF.txt file.

my %fai;

my $head;
while (<FA>) {
	s/[\r\n]+$//;
	if ($_ =~ /^>/) {
		$head = substr($_, 1);
	} else {
		$fai{$_} = $head;
#		print "$head\t$_\n";
	}
	
}

while (<TRF>)	{
	s/[\r\n]+$//;
	next if ($_ =~ /^chr/);
	my @line = split "\t", $_;
	my $double = "$line[5]$line[5]";
	my @array;
	foreach my $i (0..(length($line[5])-1)) {
		my $off = substr $double, $i, length($line[5]);
		my $RC = reverse($off);
		$RC =~ tr/ACGTacgt/TGCAtgca/;
		push @array, $off, $RC;
	}
#	print "first:\t$line[5]\n";
	foreach my $a (@array) {
#		print "$a\n";
		if ($fai{$a}) {
			print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$fai{$a}\n";
		}
	}
#	<STDIN>;
}
