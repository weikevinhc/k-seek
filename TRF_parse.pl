#!/usr/bin/env perl
# TRF_parse.pl
# parses TRF output (.dat)

use warnings;
use strict;
use POSIX;
use List::Util qw( min max );
use Math::Round;

my $dat = $ARGV[0];
my $win = 1; ##search window size. smaller size means more memory.
open DAT, "$dat" or die;
my $out = $ARGV[1] . ".txt";
my $gff = $ARGV[1] . ".gff";


my %trf;

my $chr;
while (<DAT>) {
	s/[\r\n]+$//;
	if ($_ =~ /^Seq/) {
		my @line = split " ", $_;
		$chr = $line[1];
	} elsif ($_ =~ /^\d/) {
		my @line = split " ", $_;
		my $diff = $line[1] - $line[0]  + 1;
		if ($diff >= 1000) {
			my @cor = (round($line[0]/$win), round($line[1]/$win));
			if ($trf{$chr}{$cor[0]}) {
				my $diff2 = ${$trf{$chr}{$cor[0]}}[1] -  ${$trf{$chr}{$cor[0]}}[0] + 1;
				if ($diff > $diff2) {
					foreach my $i ($cor[0] .. $cor[1]) {
					$trf{$chr}{$i} = [($line[0], $line[1], $line[1])];
					}
				}
			}
			foreach my $i ($cor[0] .. $cor[1]) {
				$trf{$chr}{$i} = [($line[0], $line[1], $line[13])];
			}
			
		}
	}
}

my @block = ("test",0 ,0);
open OUT, ">$out";
print OUT "chr\tstart\tend\tlength\tcopy#\tmotif\n";
open GFF, ">$gff";
my %trf_stats;
foreach my $i (sort keys %trf) {
	my $bcount = 0;
	my $nchar = 0;
	my $copies = 0;
	my @blockA = ();
	foreach my $c (sort {$a <=> $b} keys %{$trf{$i}}) {
		if (${$trf{$i}{$c}}[2] eq $block[0]) {
			$bcount ++;
		} else {
			if ($block[0] ne "test" && $bcount > floor(1000 / $win) ) {
				if ($bcount * $win/$nchar >= 5) {
					$trf_stats{$block[0]}{$block[1]} = $block[2];
					print OUT "$i", "\t", $blockA[0]*$win, "\t", $blockA[-1]*$win, "\t", $nchar, "\t", $bcount*$win/$nchar, "\t", "$block[0]", "\n";
					print GFF "$i\tTRF_parse\ttandemrepeats\t", $blockA[0]*$win, "\t", $blockA[-1]*$win, "\t1\t+\t.\t", "Seq $block[0]\n";
				}
			}
			$bcount = 1;
			@block = (${$trf{$i}{$c}}[2], ${$trf{$i}{$c}}[0], ${$trf{$i}{$c}}[1]);
			@blockA = ();
		}
		$nchar = length $block[0];
		$copies = floor(($block[2] - $block[1] + 1)/$nchar);
		push @blockA, $c;
		
#		print "$i\t$c\t${$trf{$i}{$c}}[0]\t${$trf{$i}{$c}}[1]\t${$trf{$i}{$c}}[2]\n";
	}
}
