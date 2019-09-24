#!/usr/bin/env perl
# TRF_blast_group.pl
# group similar tandem repeats after the self blast.

use warnings;
use strict;
use POSIX;
use List::Util qw( min max );
use List::Util qw(sum);
use Math::Round;

open FA, "$ARGV[0]" or die; ## .fa of the TRF collapsed file.
open BLAST, "$ARGV[1]" or die; ## pairwise self blast results.
## $ARGV[2] not required but should be "dup" to allow output of uncollapsed tandems

my $length = 0.9;
my $iden = 0.9;

my %fai;

my $head;
while (<FA>) {
	s/[\r\n]+$//;
	if ($_ =~ /^>/) {
		$head = substr($_, 1);
	} else {
		$fai{$head}{"length"} = length $_;
		$fai{$head}{"seq"} = $_;
#		print "$head\t$_\n";
	}
	
}

my %sim;
my %hash;
my %dup;
my %short;
my $counter = 0;
while (<BLAST>) {
	s/[\r\n]+$//;
	$counter ++;
	my @line = split "\t", $_;
	if ($fai{$line[0]}{"length"} <= 10 || $fai{$line[1]}{"length"} <= 10) {
		if ($fai{$line[0]}{"length"} <= 10) {
			$short{$line[0]} ++;
		}
		next;
	}
	$hash{$line[0]}{$line[1]}{$counter} = [@line];
}

foreach my $q (sort keys %hash) {
	foreach my $s (sort keys %{$hash{$q}}) {
		my %sub_s;
		foreach my $c (keys %{$hash{$q}{$s}}) {
			$sub_s{$c} = $hash{$q}{$s}{$c}[3];
#			print "$c = $hash{$q}{$s}{$c}[3]\n";
		}
		my @ss = sort {$sub_s{$b} <=> $sub_s{$a}} keys %sub_s;
#		print "$ss[0]\n";
		my %query;
		if (scalar @ss == 1) {
			if ($hash{$q}{$s}{$ss[0]}[3] > $length * $fai{$q}{"length"} && $hash{$q}{$s}{$ss[0]}[2] > $iden) {
				$sim{$q}{$s} ++;
			}
		} elsif (scalar @ss > 1) {
			foreach my $c (@ss) {
				next if ($hash{$q}{$s}{$c}[2] <= $iden);
#				print "seq\t";
				foreach my $i ($hash{$q}{$s}{$c}[6]..$hash{$q}{$s}{$c}[7]) {
#					print "$i";
					$query{$i} ++;
				}
#				print "\n";
			}
		}
		my %count;
		if (keys %query > $length * $fai{$q}{"length"}) {
#			print "$q\t",$fai{$q}{"seq"}, "\n", "$s\t", $fai{$s}{"seq"}, "\n";
			foreach my $k (sort {$a <=> $b} keys %query) {
#				print "$k=$query{$k}\n";
				if ($query{$k} > 1) {
					$count{"more"} ++;
				} else {
					$count{"one"} ++;
				}
			}
#			<STDIN>;
			if ($count{"more"} && $count{"more"} > 0.9 * $fai{$q}{"length"}) {
				if ($fai{$q}{"length"}*1.8 < $fai{$s}{"length"}) {
					$dup{$s}{$q} ++;
					$sim{$q}{$s} ++;
#					print "q:\t", $fai{$q}{"seq"}, "\n";
#					print "s:\t", $fai{$s}{"seq"}, "\n";
#					<STDIN>;
				}
			} elsif (keys %query > $length * $fai{$q}{"length"}) {
				$sim{$q}{$s} ++;
			}
		}
	}
}
#<STDIN>;
#print "1: ", scalar keys %sim, "\n";
#print "dup: ", scalar keys %dup, "\n";

foreach my $k (sort keys %dup) {
	print "$k\n";
	my %query;
	foreach my $c (keys %{$hash{$k}{$k}}) {
		foreach my $i ($hash{$k}{$k}{$c}[6]..$hash{$k}{$k}{$c}[7]) {
			$query{$i} ++;
		}
	}
	my %count;
	foreach my $i (sort {$a <=> $b} keys %query) {
		if ($query{$i} > 1) {
			$count{"more"} ++;
		} else {
			$count{"one"} ++;
		}
	}
	if ($count{"more"} && $count{"more"} > 0.9 *$fai{$k}{"length"}) {
		print "kept\n";
#		print "kept:\t$k\n", $fai{$k}{"seq"}, "\n";
	} else {
		delete $dup{$k};
#		print "deleted:\t$k\n", $fai{$k}{"seq"}, "\n";
	}
#	<STDIN>;
} ##reassess uncollapsed tandems
<STDIN>;
foreach my $k (sort keys %sim) {
	foreach my $s (sort keys %{$sim{$k}}) {
		unless ($sim{$s}{$k}) {
			delete $sim{$k}{$s};
		}
	}
	if ($dup{$k}) {
#		print "$k = ",$fai{$k}{"seq"}, "\n";
		delete $sim{$k};
		next;
	} else {
		foreach my $kk (sort keys %{$sim{$k}}) {
			delete $sim{$k}{kk} unless ($sim{$kk}{$k});
#			print "$kk\n";
			if ($dup{$kk}) {
				delete $sim{$k}{$kk};
			} elsif (!$sim{$kk}{$k}) {
				delete $sim{$k}{$kk};
			}
		}
	}
}
#print "2: ", scalar keys %sim, "\n";
#<STDIN>;

my %group;
my %grouped;
my $groupname = 0;
my @sim_key = sort keys %sim;
my $gcount = 0;
foreach my $k (@sim_key) {
	my %ghash = ();
	foreach my $j (keys %{$sim{$k}}) {
#		print "$k\t:\t$j\n";
		if ($grouped{$j}) {
			$ghash{$grouped{$j}} ++;
		}
	}
	if (%ghash) {
		my @sort = sort {$ghash{$b} <=> $ghash{$a}} keys %ghash;
		$grouped{$k} = $sort[0];
		$group{$sort[0]}{$k} ++;
	} else {
		$gcount ++;
#		print "$gcount\n";
		$grouped{$k} = $gcount;
		$group{$gcount}{$k} ++;
	}

}
my %gsize;
foreach my $g (sort {$a <=> $b} keys %group) {
#	print "$g\n";
	my @single = %{$group{$g}};
	$gsize{$g} = scalar keys %{$group{$g}};
	if ($gsize{$g} == 1) {
		my %ghash = ();
		foreach my $j (keys %{$sim{$single[0]}}) {
			$ghash{$grouped{$j}} ++;
		}
		if (%ghash) {
			my @sort = sort {$ghash{$b} <=> $ghash{$a}} keys %ghash;
			$grouped{$single[0]} = $sort[0];
			$group{$sort[0]}{$single[0]} ++;
			delete $group{$g};
		}
	}
}



if ($ARGV[2] && $ARGV[2] eq "dup") {
#	print "dup\n";
	foreach my $g (keys %dup) {
#		print "$g\n";
		my %ghash;
		foreach my $h (keys %{$dup{$g}}) {
			if ($grouped{$h}) {
#				print "test: $grouped{$h}\n";
				$ghash{$grouped{$h}} ++;
			}
		}
		my @sort = sort {$ghash{$b} <=> $ghash{$a}} keys %ghash;
		if ($sort[0]) {
			$group{$sort[0]}{$g} ++;
#			print "$sort[0]\n", $fai{$g}{"seq"}, "\n";
		}
#		<STDIN>;
	}
	
}

my %size;
foreach my $g (sort {$a <=> $b} keys %group) {
	my %gsize;
	foreach my $j (sort {$group{$g}{$a} cmp $group{$g}{$b}} keys %{$group{$g}}) {
		$gsize{$fai{$j}{"length"}} ++;
	}
	my @mode = sort {$gsize{$b} <=> $gsize{$a}} keys %gsize;
	$size{$g} = $mode[0];
#	print "size=$g=$mode[0]\n";
}


$counter = 0;
foreach my $g (sort {$size{$a} <=> $size{$b}} keys %group) {
	
	$counter ++;
	foreach my $j (sort {$group{$g}{$a} cmp $group{$g}{$b}} keys %{$group{$g}}) {
		print ">group$counter", "_$size{$g}", "bp_$j\n", $fai{$j}{"seq"}, "\n";
	}
}
