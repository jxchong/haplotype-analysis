#!perl
#
# Description:
#
# Usage: perl untitled.pl
#
#
# Created by Jessica on 2011-06-08

use strict;
use warnings;

if (@ARGV != 6) {
	print "Usage: $0 <phased_file> <start_bp> <stop_bp> <refhap_file> <refhap_linenum, starting from top> <pos required to be in match>\n";
	exit;
}

my $phasedfile = $ARGV[0];
my ($bpstart, $bpstop) = ($ARGV[1], $ARGV[2]);
$bpstart =~ s/,//g;
$bpstop =~ s/,//g;
my $refhapfile = $ARGV[3];
my $refhapline = $ARGV[4];
my $bpmustinclude = $ARGV[5];
$bpmustinclude =~ s/,//g;
my $maxmismatch = 3;

open (FILE, "$refhapfile") or die "Cannot open $refhapfile file.\n";
	my ($start, $stop, $mustinclude) = 0;
	my $head = <FILE>;
	$head =~ s/\s+$//;					# Remove line endings
	my @positions = split("\t", $head);
	shift(@positions);
	for (my $i=0; $i<=$#positions; $i++) {
		if ($positions[$i+1] > $bpstart) {
			$start = $i;
			last;
		}
	}
	for (my $i=0; $i<=$#positions; $i++) {
		if ($i == $#positions) {
			$stop = $i;
		} elsif ($positions[$i+1] > $bpstop) {
			$stop = $i;
			last;
		}
	}

	my @refhap;
	my $linecount = 2;
	while ( <FILE> ) {
		if ($linecount == $refhapline) {
			$_ =~ s/\s+$//;					# Remove line endings
			my @line = split ("\t", $_);	
			shift(@line);
		
			for (my $i=0; $i<=$#line; $i++) {
				if ($i >= $start && $i <= $stop) {
					push(@refhap, $line[$i]);			
				}
			}
		}
		$linecount++;
	}
close FILE;



print "subjectid\tMax match\tnSNPs checked\tmaxmatch start\tmaxmatch end\n";
open (FILE, "$phasedfile") or die "Cannot open $phasedfile file.\n";
$linecount = 0;
my ($previnfostring, $prevsubjectid, $prev_max, $prevstart, $prevstop) = ((0)x5);
<FILE>;
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	$linecount++;
	my @line = split ("\t", $_);
	my $subjectid = shift(@line);	
	my ($maxlength, $maxmatch_start, $maxmatch_end, $nzeros) = ((0)x4);
	
	for (my $offset=0; $offset<=$maxmismatch; $offset++) {
		my ($inmatch, $currlength, $currmatch_start, $currmatch_end) = ((0)x4);
		$nzeros = 0;
		my $nmismatch = $offset;
		for (my $i=0; $i<=$#line; $i++) {
			if ($i < $start) { next; }
			if ($i > $stop) {
				if ($currlength >= $maxlength) {
					$maxlength = $currlength;
					$maxmatch_start = $currmatch_start;
					$maxmatch_end = $currmatch_end;
				}
				last;
			}
	 		if ($line[$i] eq $refhap[$i] || $refhap[$i] eq '0') {
				$currlength += 1;
				if ($inmatch == 0) {
					$currmatch_start = $positions[$i];
				}
				$currmatch_end = $positions[$i];
				$inmatch = 1;	
			} elsif ($line[$i] eq '0') {
				$nzeros++;
				if ($inmatch == 1) {
					$currlength += 1;
				}
			} else {
				if ($nmismatch < $maxmismatch) {
					$nmismatch += 1;
					$currlength += 1;
					next;
				} else {
					if ($currlength >= $maxlength) {
						$maxlength = $currlength;
						$maxmatch_start = $currmatch_start;
						$maxmatch_end = $currmatch_end;
					}
					$inmatch = 0;
					$currlength = 0;
					$nmismatch = $offset;	
				}
			}
		}
	}
	
	if ($prevsubjectid == $subjectid) {
		if ($prev_max > $maxlength && ($prevstart < $bpmustinclude && $prevstop > $bpmustinclude) ) {
			print "$previnfostring\n";
		} elsif ($maxmatch_start < $bpmustinclude && $maxmatch_end > $bpmustinclude) {
			print "$subjectid\t$maxlength\t".scalar(@refhap)."\t$maxmatch_start\t$maxmatch_end\n";	
		} else {
			print "*$subjectid\tnomatch\n";
			print "*$subjectid\t$maxlength\t".scalar(@refhap)."\t$maxmatch_start\t$maxmatch_end\n";	
			print "*$previnfostring\n"; 
		}
	}
	
	
	# DEBUG
	# print "$subjectid\t$maxlength\t".scalar(@refhap)."\t$maxmatch_start\t$maxmatch_end";	
	# if ($maxmatch_start < $bpmustinclude && $maxmatch_end > $bpmustinclude) {
	# 	print "*";
	# }
	# print "\n";
	# END DEBUG
	
	$previnfostring = "$subjectid\t$maxlength\t".scalar(@refhap)."\t$maxmatch_start\t$maxmatch_end";
	$prevstart = $maxmatch_start;
	$prevstop = $maxmatch_end;
	$prevsubjectid = $subjectid;
	$prev_max = $maxlength;
}
close FILE;
