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


if (@ARGV < 4) {
	print "Usage: $0 <phased_file> <refhap_file> <refhap_linenum, starting from top> <pos of mutation> <min snps to call match> <mismatches to allow>\n";
	exit;
}






my $phasedfile = $ARGV[0];
my $refhapfile = $ARGV[1];
my $refhapline = $ARGV[2];
my $mutationbp = $ARGV[3];
$mutationbp =~ s/,//g;
my $minsnpstomatch = $ARGV[4];		# default = 100
my $maxmismatch = $ARGV[5];			# default = 2

my $countanalyzedsubj = 0;

my $mutationpos_left = 0;
my $mutationpos_right = 0;
my @refhap;


my @genotypedsubj;
open (FILE, "hutt.1415.findivlist"); 
while (<FILE>) {
	$_ =~ s/\s+$//;
	push(@genotypedsubj, $_);
}
close FILE;








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



print "subjectid\tnSNPschecked\tnMissing\tmatchedSNPs\tmaxmatch_start\tmaxmatch_end\tmatch_lengthbp\n";
open (FILE, "$phasedfile") or die "Cannot open $phasedfile file.\n";
my ($previnfostring, $prevsubjectid, $prev_max, $prev_nzeros) = ((0)x5);
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
	
	if (!grep(/^$subjectid$/, @genotypedsubj)) {
		next;
	}

	# DEBUG
	# if ($subjectid != xxx && $temp!=2) {
	# 	next;
	# } else {
	# 	$temp++;
	# 	print "$subjectid $mutationpos_left, $mutationpos_right\n";
	# }
	
	my @leftmatches;		 # store:  arrays of (nsnps, end pos of match, lengthbp, nzeros), element # is the number of mismatches
	my @rightmatches;		 # store:  arrays of (nsnps, end pos of match, lengthbp, nzeros), element # is the number of mismatches
	
	# look "before/left of" the mutation
	for (my $m=0; $m<=$maxmismatch; $m++) {
		my ($nsnps_match, $matchlengthbp, $nzeros, $nmismatch) = ((0) x 4);
		my $endmatchpos = $mutationpos_left;
		for (my $pos=$mutationpos_left; $pos>=0; $pos--) {
			if ($refhap[$pos] eq $line[$pos] || $refhap[$pos] eq '0') {
				$nsnps_match++;
				$endmatchpos = $pos;
			} elsif ($line[$pos] eq '0') {
				$nsnps_match++;
				$nzeros++;
				$endmatchpos = $pos;
			} else {
				$nmismatch++;
				$nsnps_match++;
				if ($nmismatch > $m) {
					last;
				}
			}
		}
		$leftmatches[$m] = [$nsnps_match, $positions[$endmatchpos], $positions[$endmatchpos], $nzeros];	
	}

	# look "after/right of" the mutation
	for (my $m=0; $m<=$maxmismatch; $m++) {
		my ($nsnps_match, $matchlengthbp, $nzeros, $nmismatch) = ((0) x 4);
		my $endmatchpos = $mutationpos_right;
		for (my $pos=$mutationpos_right; $pos<=$#refhap; $pos++) {
			# print "allow $m $pos $endmatchpos\n";		 # DEBUG
			if ($refhap[$pos] eq $line[$pos] || $refhap[$pos] eq '0') {
				$nsnps_match++;
				$endmatchpos = $pos;
			} elsif ($line[$pos] eq '0') {
				$nsnps_match++;
				$nzeros++;
				$endmatchpos = $pos;
			} else {
				$nmismatch++;
				$nsnps_match++;
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
	
	my $currinfostring = "$subjectid\t".scalar(@refhap)."\t$maxmatch_nzeros\t$maxsnps\t$maxmatch_start\t$maxmatch_end\t$maxlengthbp";
	
	if ($prevsubjectid == $subjectid) {
		$countanalyzedsubj++;
		# if ($maxsnps == $prev_max && $prev_max > $minsnpstomatch && ($maxmatch_nzeros/$maxsnps) <= 0.05) {
		# 	print "+=$previnfostring\n";
		# 	print "+=$currinfostring\n";
		# } els
		if (($maxsnps > $minsnpstomatch) && ($prev_max > $minsnpstomatch) && (($maxmatch_nzeros/$maxsnps) <= 0.05) && (($prev_nzeros/$prev_max) <= 0.05)) {
			print "++$previnfostring\n";
			print "++$currinfostring\n";
		} elsif ($prev_max > $maxsnps && $prev_max > $minsnpstomatch) {
			if ($prev_nzeros/$prev_max <= 0.05) {
				print "$previnfostring\n";
			} else {
				print "?$previnfostring\n";
			}
		} elsif ($maxsnps > $prev_max && $maxsnps > $minsnpstomatch) {
			if ($maxmatch_nzeros/$maxsnps <= 0.05) {
				print "$currinfostring\n";
			} else {
				print "?$currinfostring\n";
			}
		} else {
			print "*$subjectid\tnomatch\n";
			print "*$previnfostring\n";
			print "*$currinfostring\n";
		}
	}


	# DEBUG 
	# print "$subjectid\t".scalar(@refhap)."\t$maxsnps\t$maxmatch_start\t$maxmatch_end\t$maxlengthbp";	
	# print "\n";

	
	$previnfostring = "$subjectid\t".scalar(@refhap)."\t$maxmatch_nzeros\t$maxsnps\t$maxmatch_start\t$maxmatch_end\t$maxlengthbp";
	$prevsubjectid = $subjectid;
	$prev_max = $maxsnps;
	$prev_nzeros = $maxmatch_nzeros;
}
close FILE;

print STDERR "Analyzed $countanalyzedsubj subjects with genotype information\n";
