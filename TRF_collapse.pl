#!usr/bin/perl
# TRF_collapse.pl
# collapse trf output into a fastA

use warnings;
use strict;
use POSIX;
use List::Util qw( min max );
use Math::Round;

my $dat = $ARGV[0];
open DAT, "$dat" or die;

my %trf;


my $chr;
my %hash;
while (<DAT>) {
	s/[\r\n]+$//;
	if ($_ =~ /^Seq/) {
		my @line = split " ", $_;
		$chr = $line[1];
	} elsif ($_ =~ /^\d/) {
		my @line = split " ", $_;
		my $diff = $line[1] - $line[0]  + 1;
		if ($diff >= 1000 && $line[3] >= 10) {
			my $kmer = &RC_sort($line[13]);
			$hash{$kmer} = length($kmer);
		}
	}
}

my $count = 0;
foreach my $k (sort {$hash{$a} <=> $hash{$b}} keys %hash) {
	$count ++;
	print ">TRF$count","_$hash{$k}\n";
	print "$k\n";
}

sub RC_sort {
	my $double = "$_[0]$_[0]";
	my @array;
	foreach my $i (0..(length($_[0])-1)) {
		my $off = substr $double, $i, length($_[0]);
		my $RC = reverse($off);
		$RC =~ tr/ACGTacgt/TGCAtgca/;
		push @array, $off, $RC;
	}
	@array = sort @array;
	return $array[0];
}