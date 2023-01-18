#!usr/bin/perl
# k_seek.r4.pl
# identifies repeats in every read, and counts them, version 4 (combines k_finder and k_counter).
use strict;
use warnings;
use POSIX;

my $input1 = shift @ARGV;
my $input2 = shift @ARGV;
open GENOME_DATA, $input1 or die;
open OUTPUT , ">$input2.rep";
open TOTAL, ">$input2.rep.total";

my $seq_ID = "undeclared";
my $seq = "undeclared";
my $seq_line3 = "undeclared";
my $quality = "undeclared";

my @prime = (1,2,3,5,7,11);
my %cutoff = (
	"1" => 0, "2" => 5, "3" => 3,
	"4" => 2, "5" => 2, "6" => 2,
	"7" => 1, "8" => 1, "9" => 1, 
	"10" => 1, "11" => 1, "12" => 1, 
	"13" => 1, "14" => 1, "15" => 1,
	"16" => 1, "17" => 1, "18" => 1,
	"19" => 1, "20" => 1
); # cut-off for the minimium number of values for each kmer stored in hash

my $line_counter = 0;
my %kmer_total = ();
READS: while (<GENOME_DATA>) {
	s/[\r\n]+$//;
	$line_counter ++;
	if ($line_counter == 1) {
		$seq_ID = $_ ;
	} elsif ($line_counter == 2) {
		$seq = $_ ;
	} elsif ($line_counter == 3) {
		$seq_line3 = $_ ;
	} elsif ($line_counter == 4) {
		my $hit = 0;
		$line_counter = 0;
		$quality = $_;
		
		
		if ((length $seq < 20) or (substr($seq, -5, 5) eq "NNNNN")) {
#			print length($seq), "\n";
#			print substr($seq, -5, 5), "\n";
#			<STDIN>;
			next;
		}
#		print "$seq_ID\n$seq\n$seq_line3\n$quality\n";
		
		my %lengthHoH;
		my $href;
		SEQ: foreach my $k (reverse(6..11)) {
			$href = { rep_identify($seq, $k) };
			next if (keys %{$href} < 1 or keys %{$href} > 3);
			foreach my $key (keys %{$href}) {
#				print "$k mer: $key with ${$href}{$key} repeats\n";
				if ($k * ${$href}{$key} >= 25) {
					my ($sub1, $sub2) = internal($key, ${$href}{$key});
					if ($sub1) {
						delete ${$href}{$key};
						${$href}{$sub1} = $sub2;
#						print "internal: $sub1 = $sub2\n";						
					}
					$lengthHoH{$k} = $href;
				}
			}
		} ## put all quantified length into hoh
		if (keys %lengthHoH >= 1) {
#			foreach my $k (sort { $a <=> $b} keys %lengthHoH) {
#				print scalar keys %{$lengthHoH{$k}}, ":\t";			
#				foreach my $key ( keys %{$lengthHoH{$k}} ) {
#					print "$key = ${$lengthHoH{$k}}{$key}\t";
#				}
#				print "\n";			
#			}
			
			my %sizeHoH;
			
			foreach my $k (keys %lengthHoH) {
				foreach my $key (keys %{$lengthHoH{$k}}) {
#					print "$k\t$key =", ${$lengthHoH{$k}}{$key} * length($key),"\n";
					if ($sizeHoH{$key}) {
						if ($sizeHoH{$key} < ${$lengthHoH{$k}}{$key} * length($key)) {
							$sizeHoH{$key} = ${$lengthHoH{$k}}{$key} * length($key);
						}
					} else {
						$sizeHoH{$key} = ${$lengthHoH{$k}}{$key} * length($key);
					}					
				}
			}
			my @sorted = sort {$sizeHoH{$b} <=> $sizeHoH{$a}} keys %sizeHoH;
			
#			foreach my $k (@sorted) {
#				print "$k = $sizeHoH{$k}\n";
#			}

			my $kmer = $sorted[0];
			my ($kmernew, $kcount) = kcount($kmer, $seq);
			$kmer = $kmernew;
			if ($kcount * length($kmer) > 50) {
				
				print OUTPUT "$seq_ID\n$seq\n$seq_line3\n$quality\n";
				print OUTPUT "$kmer=$kcount\n";
				$kmer_total{$kmer} += $kcount; # logs and adds the repeat to the total hash
				$hit += 1;
				
			}
		}
#		print "$seq_ID\n$seq\n";
#		print "$sorted[0]=",$sizeHoH{$sorted[0]}/length($sorted[0]),"\n\n";
#		<STDIN>;
		
		unless ($hit) {
			my $kmer;
			my %lengthHoH_long;
			my $href_long;
			foreach my $k (reverse(12..20)) {
				$href_long = { rep_identify($seq, $k) };
				if (keys %{$href_long} >= 1 or  %{$href_long} <= 2) {
					foreach my $key (keys %{$href_long}) {
	#					print "$k mer: $key with ${$href_long}{$key} repeats\n";
						if ($k * ${$href_long}{$key} >= 30) {
							my $keykey = "$key$key";
							if ($seq =~ /$keykey/) {
								my ($sub1, $sub2) = internal($key, ${$href_long}{$key});
							 	if ($sub1) {
#									print "$key = ${$href_long}{$key} => $sub1 = $sub2\n";
									delete ${$href_long}{$key};
								} else {	
									$lengthHoH_long{$k} = $href_long;
								}
							}
						}
					}
				}	
			}
			if (keys %lengthHoH_long == 1 or keys %lengthHoH_long == 2) {
				my %hash = ();
				foreach my $k (keys %lengthHoH_long) {
					if (keys %{$lengthHoH_long{$k}} == 1 or keys %{$lengthHoH_long{$k}} == 2) {						
						foreach my $l (sort { $lengthHoH_long{$k}{$b} <=> $lengthHoH_long{$k}{$a} } keys %{$lengthHoH_long{$k}}) {
							$hash{$l} = $lengthHoH_long{$k}{$l} * $k;
						}
					}
				}
				my @sorted = (sort { $hash{$b} <=> $hash{$a} } keys %hash);
				$kmer = $sorted[0];
			}
			if ($kmer) {
				my ($kmernew, $kcount) = kcount($kmer, $seq);
				$kmer = $kmernew;
				if ($kcount * length($kmer) > 50) {				
					print OUTPUT "$seq_ID\n$seq\n$seq_line3\n$quality\n";
					print OUTPUT "$kmer=$kcount\n"; 
					$kmer_total{$kmer} += $kcount; # logs and adds the repeat to the total hash 
				}
			}
		}
	}
}

foreach my $rep (sort keys %kmer_total) {
	print TOTAL "$rep\t$kmer_total{$rep}\n";
}



sub rep_identify {
	## kmer length is $_[0]
	my $string = $_[0];
	my $length = $_[1];
	my @rep_array = ($string =~ m/(.{1,$length})/gs); 
	
	my %rep_hash = ();
	for (@rep_array) {
#		print "$_\t";
		$rep_hash{$_} ++;
	}
	foreach my $kmer (keys %rep_hash) {
		delete $rep_hash{$kmer} if ($rep_hash{$kmer} <= $cutoff{$length});
		
	}
	
	if (keys %rep_hash == 1) { # if only one element in hash, then it is logged as repeat
		foreach my $kmer (keys %rep_hash) {
			return %rep_hash;
		}
	} elsif (keys %rep_hash > 1 && keys %rep_hash <= 3) { # if more than one element in hash, check for offsets
		my @karray = sort {$rep_hash{$a} <=> $rep_hash{$b}} keys %rep_hash;
		my $test = shift @karray;
#		print "test: $test : $rep_hash{$test}\n";
		foreach my $kmer (@karray) {
#			print "multi: $kmer : $rep_hash{$kmer}\n";
			if (testoff($test, $kmer)) {
				$rep_hash{$test} += $rep_hash{$kmer};
#				print "delete: $kmer ($test): $rep_hash{$kmer}\n";
				delete $rep_hash{$kmer};
			}
		}
		foreach my $kmer (keys %rep_hash) {
#			print "remaining: $kmer : $rep_hash{$kmer}\n";
		}
		return %rep_hash;
	} else {
		return %rep_hash;
	}
}



sub internal { 
	# check for internal repetition in the identified repeat
	## repeat to be tested is $_[0], it's occurence is $_[1];
	
		
	my $rep = $_[0];
	my $count = $_[1];
	my %return_hash;
	$return_hash{$rep} = $count;
	my %mono_hash;
	my @mono_array = ($rep =~ m/(.{1,1})/gs);
	foreach (@mono_array) {
		$mono_hash{$_} ++;
	}
	if (keys %mono_hash == 1) {
		return($mono_array[0], $count*length($rep));
	} elsif (grep { $_ == length($rep) } @prime) {
		return(0,0);
	} else {
		my @factors = grep { length($rep) % $_ == 0 } 2 .. (floor(length($rep))/2);
		INTER: foreach my $factor (@factors) {
			my @factor_array = ($rep =~ m/(.{1,$factor})/gs); 
			my %factor_hash;
			foreach (@factor_array) {
				$factor_hash{$_} ++;
			}
			if ( keys %factor_hash == 1 ) {
				foreach my $key (keys %factor_hash) {
					$factor_hash{$key} = $count*length($rep)/$factor;
#					print "$seq_ID\n$seq\n$seq_line3\n$quality\n";
#					print "hit!\t $rep = $count\t$key = ", $count*length($rep)/$factor, "\n";
#					<STDIN>;
					return($key, $count*length($rep)/$factor);					
				}
			}
		}
	}
	return (0,0);
}

sub testoff
{

	if ($_[0] eq $_[1]) {
		return(1);
	} else {
		my $r2 = "$_[0]$_[0]";
		if ($r2 =~ /$_[1]/) {
#			print "$_[0]\t$_[1]\n";
			return(1);
		} else {
			return(0);
		}
	}
}

sub degen { 
	# identify single substitution/insertion/deletion within a kmer
	my $target = $_[0];
	my $kmer2 = $_[1];
	my $diff = length($target) - length($kmer2);
	if ($diff == 0) {
		# find substitutions
		# kmer length is positive doing 1 - bigger than one yield a negative number
		foreach (0 .. ( length($kmer2) - 1)) {
			# unfortunalty the subtr function return the exised string not the new string and in fact modify the original string
			my $temp_k = $kmer2; # make a copy copy
			my $degen_kmer = substr($temp_k, $_, 1, "[ACTGN]"); 
			if ($target =~ /$temp_k/g) {
				return 1;
			}
		} 
		return 0;
	} elsif ($diff > 0) {
		# find insertions
		foreach (0..length($kmer2)) {
			my $tar_kmer = substr($target, 0, $_) . substr($target, $_ + $diff);
			if ($tar_kmer eq $kmer2) {
				return 1;
			}
		}
		return 0;
	} elsif ($diff < 0 && length($target) >= 3) {
		# find deletions
		foreach (0..length($target)) {
			my $del_kmer = substr($kmer2, 0, $_) . substr($kmer2, $_ - $diff);
			if ($target eq $del_kmer) {
				return 1;
			}
		}
		return 0;
	} else {
		return 0;
	}			
}

sub kcount {
	my @sub_seq = ();
	my $k = $_[0];
	my $seq_copy = $_[1];
	my $string_counter = 0;
	if (length($k) <= 3) {
		$string_counter = index($seq_copy, "$k$k");
	} else {
		$string_counter = index($seq_copy, "$k");
	}
	if ($string_counter) {
		foreach my $i (1..(length($k)-1)) {
			if (substr($seq_copy, $string_counter - 1, 1) eq substr($k, -1)) {
				$k = substr($seq_copy, $string_counter - 1, length($k));
				$string_counter -= 1;
			} else {
				last;
			}
		}
	} #check for the best offset of the kmer.
	
	while ($string_counter != -1) {
		if ($string_counter > 0) {
			push @sub_seq, substr($seq_copy, 0, $string_counter);
			$seq_copy = substr($seq_copy, $string_counter);
		} elsif ($string_counter == 0) {
			push @sub_seq, substr($seq_copy, 0, length($k));
			$seq_copy = substr($seq_copy, length($k));
		}
		$string_counter = index($seq_copy, $k);
	}
	if ($seq_copy) {
		push @sub_seq, $seq_copy;
	}
	unshift @sub_seq, "@@";
	push @sub_seq, "%%"; ## adding random sequences at front and back of array
	
	my $kcounts = 0;
	my $degcounts = 0;
	foreach my $i (1..(scalar(@sub_seq)-2)) {
		if ($sub_seq[$i] eq $k) {
			if ($sub_seq[$i] eq $sub_seq[$i-1] or $sub_seq[$i] eq $sub_seq[$i+1]) {
				$kcounts ++;
			}
		} elsif ($sub_seq[$i] ne $k) {
			if ($sub_seq[$i-1] eq $sub_seq[$i+1]) {
				if (degen($sub_seq[$i], $k)) {
					$kcounts ++;
					$degcounts ++;
				}
					
			}
		}
	}
		
	return($k, $kcounts);
}
	
