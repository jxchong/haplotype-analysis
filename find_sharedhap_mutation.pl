#!perl
#
# Description: Find longest shared homozygous haplotype (vs reference haplotype) surrounding mutation position based on IBS>=1 sharing
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


my @badfindivs;
open (FILE, "MendErr.2012-01-17.findivlist"); 
while (<FILE>) {
	$_ =~ s/\s+$//;
	unless ($_ == 173142 || $_ == 171351 || $_ == 173152) {
		push(@badfindivs, $_);
		
	}
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
my $mutationpos_gap_bp = $positions[$mutationpos_right] - $positions[$mutationpos_left];


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





print "subjectid\tnSNPschecked\tnMissing\tmatchedSNPs\tmaxmatch_start_snpno\tmaxmatch_end_snpno\tmaxmatch_start_bp\tmaxmatch_end_bp\tmatch_lengthbp\n";
open (FILE, "$phasedfile") or die "Cannot open $phasedfile file.\n";
my (@previnfo, @prevhaplotype);
my ($prevsubjectid, $prev_max, $prev_nzeros) = ((0)x5);
my $consensus_start = my $consensus_end = 0;			# start SNP position number, end SNP position number
my $temp = 1;
<FILE>;
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @currenthaplotype = split ("\t", $_);
	my $subjectid = shift(@currenthaplotype);

	if (scalar(@currenthaplotype) != scalar(@refhap)) {
		print STDERR "Error, refhap and phased haps have different number of SNPs\n";
		die;
	}
	
	if (!grep(/^$subjectid$/, @genotypedsubj) || (grep(/^$subjectid$/, @badfindivs) && $phasedfile !~ 'SMN')) {
		next;
	}

	if ($prevsubjectid == $subjectid) {
		my @leftmatches;		 # store:  arrays of (nsnps, end snp pos of match, lengthbp, nzeros), element # is the number of mismatches
		my @rightmatches;		 # store:  arrays of (nsnps, end snp pos of match, lengthbp, nzeros), element # is the number of mismatches

		# look "before/left of" the mutation
		for (my $m=0; $m<=$maxmismatch; $m++) {
			my ($nsnps_match, $matchlengthbp, $nzeros, $nmismatch) = ((0) x 4);
			my $endmatchpos = $mutationpos_left;
			for (my $pos=$mutationpos_left; $pos>=0; $pos--) {
				if (($refhap[$pos] eq $currenthaplotype[$pos] && $refhap[$pos] eq $prevhaplotype[$pos]) || $refhap[$pos] eq '0') {
					$nsnps_match++;
					$endmatchpos = $pos;
				} elsif ($currenthaplotype[$pos] eq '0' || $prevhaplotype[$pos] eq '0') {
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
			$leftmatches[$m] = [$nsnps_match, $endmatchpos, ($positions[$endmatchpos]-$positions[$mutationpos_left]), $nzeros];	
		}

		# look "after/right of" the mutation
		for (my $m=0; $m<=$maxmismatch; $m++) {
			my ($nsnps_match, $matchlengthbp, $nzeros, $nmismatch) = ((0) x 4);
			my $endmatchpos = $mutationpos_right;
			for (my $pos=$mutationpos_right; $pos<=$#refhap; $pos++) {
				# print "allow $m $pos $endmatchpos\n";		 # DEBUG
				if (($refhap[$pos] eq $currenthaplotype[$pos] && $refhap[$pos] eq $prevhaplotype[$pos]) || $refhap[$pos] eq '0') {
					$nsnps_match++;
					$endmatchpos = $pos;
				} elsif ($currenthaplotype[$pos] eq '0' || $prevhaplotype[$pos] eq '0') {
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
			$rightmatches[$m] = [$nsnps_match, $endmatchpos, ($positions[$endmatchpos]-$positions[$mutationpos_right]), $nzeros];
		}

		# check all possible combinations of matching segments to find the longest (based on bp length)
		my ($maxsnps, $maxlengthbp, $maxmatch_start_snpno, $maxmatch_end_snpno, $maxmatch_nzeros, $maxmatch_nmismatch, $maxmatch_start_bp, $maxmatch_end_bp) = ((0)x8);

		for (my $l=0; $l<=$maxmismatch; $l++) {
			my $r = ($maxmismatch-$l);
			my @leftmatch = @{$leftmatches[$l]};		# array of (nsnps, end pos of match, lengthbp, nzeros)
			my @rightmatch = @{$rightmatches[$r]};	

			my $currlengthbp = $rightmatch[2]-$leftmatch[2]+$mutationpos_gap_bp;
			my $currsnpsmatch;
			if ($mutationpos_left != $mutationpos_right) {
				$currsnpsmatch = $leftmatch[0]+$rightmatch[0];
			} else {
				$currsnpsmatch = $leftmatch[0]+$rightmatch[0]-1;
			}

			if ($currsnpsmatch > $maxsnps) {
				$maxlengthbp = $currlengthbp;
				$maxsnps = $currsnpsmatch;
				$maxmatch_start_snpno = $leftmatch[1];
				$maxmatch_end_snpno = $rightmatch[1];
				$maxmatch_start_bp = $positions[$leftmatch[1]];
				$maxmatch_end_bp = $positions[$rightmatch[1]];
				$maxmatch_nzeros = $leftmatch[3]+$rightmatch[3];
			}			
		}
				
		my @currinfo = ($subjectid, scalar(@refhap), $maxmatch_nzeros, $maxsnps, $maxmatch_start_snpno, $maxmatch_end_snpno, $maxmatch_start_bp, $maxmatch_end_bp, $maxlengthbp);
		
		$countanalyzedsubj++;
		if (($maxsnps > $minsnpstomatch) && (($maxmatch_nzeros/$maxsnps) <= 0.05)) {
			print join("\t", @currinfo)."\n";
		} elsif ($maxsnps > $minsnpstomatch) {
			print "?".join("\t", @currinfo)."\n";
		} else {
			# print "*$subjectid\tnomatch\n";
			print "*".join("\t", @currinfo)."\n";
		}
	}
	
	$prevsubjectid = $subjectid;
	# @previnfo = ($subjectid, scalar(@refhap), $maxmatch_nzeros, $maxsnps, $maxmatch_start_snpno, $maxmatch_end_snpno, $maxmatch_start_bp, $maxmatch_end_bp, $maxlengthbp);
	# $prev_max = $maxsnps;
	# $prev_nzeros = $maxmatch_nzeros;
	@prevhaplotype = @currenthaplotype;
}
close FILE;

print STDERR "Analyzed $countanalyzedsubj subjects with genotype information, allowing $maxmismatch mismatches around mutation (SNPno $mutationpos_left $mutationpos_right)\n";
