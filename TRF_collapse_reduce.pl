#!usr/bin/perl
# TRF_collapse_reduce.pl
# parses TRF output (.dat)

use warnings;
use strict;
use POSIX;
use List::Util qw( min max );
use Math::Round;

my $dat = $ARGV[0];
open DAT, "$dat" or die;
open OUT1, ">$ARGV[1].recol.fa";
open OUT2, ">$ARGV[1].uncol.fa";
my %trf;
my %seq;
my $chr;
while (<DAT>) {
	s/[\r\n]+$//;
	if ($_ =~ /^Seq/) {
		my @line = split " ", $_;
		$chr = $line[1];
	} elsif ($_ =~ /^\d/) {
		my @line = split " ", $_;
		next if ($line[1] - $line[0] < 20);
		my $kmer = &RC_sort($line[13]);
		$trf{$chr}{$kmer} = $line[1] - $line[0];
		$seq{$chr}{$kmer} = $line[14];
	}
}

my %trf_all;
my $count = 0;
foreach my $chr (keys %trf) {
	my @k = sort {$trf{$chr}{$b} <=> $trf{$chr}{$a}} keys %{$trf{$chr}};
	unless ($trf_all{$k[0]}) {
		$count ++;
		print OUT1 ">TRF_collapsed$count\n";
		print OUT1 "$k[0]\n";
#		print "$k[0]\n";
#		print "$seq{$chr}{$k[0]}\n";
		if (length($k[0])*1.8 < length($seq{$chr}{$k[0]}) ) {
			print OUT2 ">TRF_collapsed$count", "_", length($seq{$chr}{$k[0]}), "\n";
			print OUT2 "$seq{$chr}{$k[0]}\n";
		}
		$trf_all{$k[0]} ++;
	}
	
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