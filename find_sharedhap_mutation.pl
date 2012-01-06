#!perl
#
# Description: Find longest shared haplotype (vs reference haplotype) surrounding mutation position.
#
# Usage: perl untitled.pl
#
#
# Created by Jessica on 2012-01-04

use strict;
use warnings;


if (@ARGV != 4) {
	print "Usage: $0 <phased_file> <refhap_file> <refhap_linenum, starting from top> <pos required to be in match>\n";
	exit;
}

my $phasedfile = $ARGV[0];
my $refhapfile = $ARGV[1];
my $refhapline = $ARGV[2];
my $mutationbp = $ARGV[3];
$mutationbp =~ s/,//g;
my $maxmismatch = 3;


my $mutationpos_left = 0;
my $mutationpos_right = 0;
my @refhap;

open (FILE, "$refhapfile") or die "Cannot open $refhapfile file.\n";
my $head = <FILE>;
$head =~ s/\s+$//;					# Remove line endings
my @positions = split("\t", $head);
shift(@positions);
for (my $i=0; $i<=$#positions; $i++) {
	if ($positions[$i] < $mutationbp && $positions[$i+1] > $mutationbp) {
		$mutationpos_left = $i;
		$mutationpos_right = $i+1;
	}
	if ($positions[$i] == $mutationbp) {
		$mutationpos_left = $mutationpos_right = $i;
	}
}


my $linecount = 2;
while ( <FILE> ) {
	if ($linecount == $refhapline) {
		$_ =~ s/\s+$//;					# Remove line endings
		my @line = split ("\t", $_);	
		shift(@line);
		@refhap = @line;
		last;
	}
	$linecount++;
}
close FILE;



print "subjectid\tnSNPschecked\tmatchedSNPs\tmaxmatch_start\tmaxmatch_end\tmatch_lengthbp\n";
open (FILE, "$phasedfile") or die "Cannot open $phasedfile file.\n";
my ($previnfostring, $prevsubjectid, $prev_max) = ((0)x5);
my $temp = 1;
<FILE>;
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	my $subjectid = shift(@line);

	
	if (scalar(@line) != scalar(@refhap)) {
		print "Error, refhap and phased haps have different number of SNPs\n";
		die;
	}
	
	my @leftmatches;		 # store:  arrays of (nsnps, end pos of match, lengthbp, nzeros), element # is the number of mismatches
	my @rightmatches;		 # store:  arrays of (nsnps, end pos of match, lengthbp, nzeros), element # is the number of mismatches
	
	# look "before/left of" the mutation
	for (my $m=0; $m<=$maxmismatch; $m++) {
		my ($nsnps_match, $endmatchpos, $matchlengthbp, $nzeros, $nmismatch) = ((0) x 5);
		for (my $pos=$mutationpos_left; $pos>=0; $pos--) {
			if ($refhap[$pos] eq $line[$pos] || $refhap[$pos] eq '0') {
				# print "at $pos $positions[$pos], match of $refhap[$pos] with $line[$pos]\n";
				$nsnps_match++;
				$endmatchpos = $pos;
			} elsif ($line[$pos] eq '0') {
				$nsnps_match++;
				$nzeros++;
				$endmatchpos = $pos;
			} else {
				$nmismatch++;
				# print "at $pos $positions[$pos], no match of $refhap[$pos] with $line[$pos].   $nsnps_match, $positions[$endmatchpos], $positions[$endmatchpos], $nzeros\n";
				if ($nmismatch > $m) {
					last;
				}
			}
		}
		$leftmatches[$m] = [$nsnps_match, $positions[$endmatchpos], $positions[$endmatchpos], $nzeros];	
	}

	# look "after/right of" the mutation
	for (my $m=0; $m<=$maxmismatch; $m++) {
		my ($nsnps_match, $endmatchpos, $matchlengthbp, $nzeros, $nmismatch) = ((0) x 5);
		for (my $pos=$mutationpos_right; $pos<=$#refhap; $pos++) {
			if ($refhap[$pos] eq $line[$pos] || $refhap[$pos] eq '0') {
				$nsnps_match++;
				$endmatchpos = $pos;
			} elsif ($line[$pos] eq '0') {
				$nsnps_match++;
				$nzeros++;
				$endmatchpos = $pos;
			} else {
				$nmismatch++;
				if ($nmismatch > $m) {
					last;
				}
			}
		}
		$rightmatches[$m] = [$nsnps_match, $positions[$endmatchpos], $positions[$endmatchpos], $nzeros];
	}

	# check all possible combinations of matching segments to find the longest (based on bp length)
	my ($maxsnps, $maxlengthbp, $maxmatch_start, $maxmatch_end, $maxmatch_nzeros) = ((0)x5);
		
	for (my $l=0; $l<=$maxmismatch; $l++) {
		my $r = ($maxmismatch-$l);
		my @leftmatch = @{$leftmatches[$l]};		# array of (nsnps, end pos of match, lengthbp, nzeros)
		my @rightmatch = @{$rightmatches[$r]};	
		
		my $currlengthbp = $rightmatch[2]-$leftmatch[2];
		my $currsnpsmatch;
		if ($mutationpos_left != $mutationpos_right) {
			$currsnpsmatch = $leftmatch[0]+$rightmatch[0];
		} else {
			$currsnpsmatch = $leftmatch[0]+$rightmatch[0]-1;
		}
	
		if ($currsnpsmatch > $maxsnps) {
			$maxlengthbp = $currlengthbp;
			$maxsnps = $currsnpsmatch;
			$maxmatch_start = $leftmatch[1];
			$maxmatch_end = $rightmatch[1];
			$maxmatch_nzeros = $leftmatch[3]+$rightmatch[3];
		}			
	}
	

	print "$subjectid\t".scalar(@refhap)."\t$maxsnps\t$maxmatch_start\t$maxmatch_end\t$maxlengthbp";	
	print "\n";

	
	$previnfostring = "$subjectid\t".scalar(@refhap)."\t$maxsnps\t$maxmatch_start\t$maxmatch_end\t$maxlengthbp";
	$prevsubjectid = $subjectid;
	$prev_max = $maxlengthbp;
}
close FILE;