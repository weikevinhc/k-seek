#!usr/bin/perl
# k_compiler.pl
# colapses all rep.total files in directory into one table
use warnings;
use strict;
my $input = shift @ARGV;
my $input2 = shift @ARGV;
opendir (DIR, $input) or die $!;

my %rep_hash = ();
my @files;
while (my $file = readdir(DIR)) {
	next if ($file =~ m/^\./);
	print "$input\/$file\n";
	push @files, $file;
	open REP_TOTAL, "$input\/$file" or die;
	while (<REP_TOTAL>) {
		s/[\r\n]+$//;
		my @split = split "\t", $_;
		my $kmer = $split[0];
		my $kmerx2 = "$kmer$kmer";
		my @array = ();
		my $sort_counter = 0;
		my $kmer_length = length($kmer);
		while ($sort_counter < $kmer_length) {
			my $subkmer = (substr $kmerx2, $sort_counter, $kmer_length);
			push @array, $subkmer;
			my $revcomp = reverse($subkmer);
			$revcomp =~ tr/ACGTacgt/TGCAtgca/;
			push @array, $revcomp;
			$sort_counter ++;
		}
		@array = sort(@array); # alphabetizing all possible offsets of the repeat.
		$kmer = $array[0];
		if ($rep_hash{$kmer}) {
			next;
		} else {
			$rep_hash{$kmer} = 0;
		}
	}
	close REP_TOTAL;
}
closedir (DIR);

opendir (DIR, $input) or die $!;

my @rep_array;
foreach my $rep (keys %rep_hash) {
	push @rep_array, $rep;
}
@rep_array = sort @rep_array;
open OUTPUT, ">$input2.rep.compiled";
print OUTPUT "lines";
my @rep_length_array;
foreach (@rep_array) {
	my $revcomp = reverse($_);
	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
	print OUTPUT "\t", $_, "/", $revcomp;
	push @rep_length_array, length $_;	
}
print OUTPUT "\ttotal_bp\n";

foreach my $file (sort @files) {
	print "$input\/$file counting\n";
	foreach my $rep (keys %rep_hash) {
		$rep_hash{$rep} = 0;
	}
	open REP_TOTAL, "$input\/$file" or die;
	while (<REP_TOTAL>) {
		s/[\r\n]+$//;
		my @split = split "\t", $_;
		my $kmer = $split[0];
		my $kmerx2 = "$kmer$kmer";
		my @array = ();
		my $sort_counter = 0;
		my $kmer_length = length($kmer);
		while ($sort_counter < $kmer_length) {
			my $subkmer = (substr $kmerx2, $sort_counter, $kmer_length);
			push @array, $subkmer;
			my $revcomp = reverse($subkmer);
			$revcomp =~ tr/ACGTacgt/TGCAtgca/;
			push @array, $revcomp;
			$sort_counter ++;
		}
		@array = sort(@array);
		$kmer = $array[0];
		$rep_hash{$kmer} = $rep_hash{$kmer} + $split[1];
	}
	
	print OUTPUT "$file\t";
	my $total = 0;
	foreach my $kmer (@rep_array) {
		my $kmerx2 = "$kmer$kmer";
		my @array = ();
		my $sort_counter = 0;
		my $kmer_length = length($kmer);
		while ($sort_counter < $kmer_length) {
			my $subkmer = (substr $kmerx2, $sort_counter, $kmer_length);
			push @array, $subkmer;
			my $revcomp = reverse($subkmer);
			$revcomp =~ tr/ACGTacgt/TGCAtgca/;
			push @array, $revcomp;
			$sort_counter ++;
		}
		my $kmer_total = 0;
		foreach my $rep (@array) {
			if ($rep_hash{$rep}) {
				$kmer_total += $rep_hash{$rep};
				$rep_hash{$rep} = 0;
			} else {
				next;
			}
		}
		print OUTPUT $kmer_total, "\t";
		my $rep_length = length $kmer;
		my $rep_bases = $kmer_total * $rep_length;
		$total += $rep_bases;
	}
	print OUTPUT $total, "\n";
	close REP_TOTAL;
}
		