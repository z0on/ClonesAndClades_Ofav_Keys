#!/usr/bin/env perl

my $usage="

zooxType2.pl:

Calculates relative representations of zooxanthellae clades (A,B,C,D)
using SAM files with RAD reads mapped to combined reference in which
host: chr1-chr10 (only matches to chromosome 1 will be counted)
A: chr11
B: chr12
C: chr13
D: chr14

Calculates number of high-quality (= high uniqueness) mappings a proxy of relative 
clade abundance in each SAM file. 

arguments and defaults:

    ext=sam : extension of SAM files
    minq=40 : minimum mapping quality
host=\"chr1\" : name of host contig to use for calculating zoox prevalence

output: 

tab-delimited table (printed to STDOUT) of relative proportions of clades,
in the form of 

sample\tsitesHost\tsitesA\tsitesB\tsitesC\tsitesD\tcountHost\tcountA\tcountB\tcountC\tcountD\tfracZoox\tfracA\tfracB\tfracC\tfracD

";

my $ext="sam";
if (" @ARGV "=~/glob=(\S+)/) { $ext=$1;}
my $minq=40;
if (" @ARGV "=~/minq=(\d+)/) { $minq=$1;}
my $host="chr1";
if (" @ARGV "=~/host=(\S+)/) { $host=$1;}

opendir THIS, ".";
my @files=grep/\.$ext$/, readdir THIS;
@files=sort {$a cmp $b} @files;
if ($#files==-1) { die $usage; }

print "sample\tsitesH\tsitesA\tsitesB\tsitesC\tsitesD\tcountH\tcountA\tcountB\tcountC\tcountD\tfracZ\tfracA\tfracB\tfracC\tfracD\n";

foreach my $file (@files){
	my $sH=0;
	my $sA=0;
	my $sB=0;
	my $sC=0;
	my $sD=0;
	my $nH=0;
	my $nA=0;
	my $nB=0;
	my $nC=0;
	my $nD=0;
	my $fZ=0;
	my $fA=0;
	my $fB=0;
	my $fC=0;
	my $fD=0;
	my %seenH={};
	my %seenA={};
	my %seenB={};
	my %seenC={};
	my %seenD={};
	my $site;
	my $qual;
	next if ($file=~/HASH/);
	open INF, $file or die "$usage\n\nCannot open sam file $file\n\n";
	
	while (<INF>) {
                if ($_=~/\t$host\t(\d+)\t(\d+)/) {
                        $qual=$2;
                        $site=$1;
                        next if ($qual<$minq);
#warn "$file:A:Q=$qual:S=$site\n$_";
                        $nH++;
                        if (!$seenH{$site}) {
                                $seenH{$site}=1;
                                $sH++;
                        }
                }
		elsif ($_=~/\tchr11\t(\d+)\t(\d+)/) {
			$qual=$2;
			$site=$1;
			next if ($qual<$minq);
#warn "$file:A:Q=$qual:S=$site\n$_"; 
			$nA++;
			if (!$seenA{$site}) {
				$seenA{$site}=1;
				$sA++;
			}
		}
		elsif ($_=~/\tchr12\t(\d+)\t(\d+)/) {
			$qual=$2;
			$site=$1;
			next if ($qual<$minq);
			$nB++;
			if (!$seenB{$site}) {
				$seenB{$site}=1;
				$sB++;
			}
		}
		elsif ($_=~/\tchr13\t(\d+)\t(\d+)/) {
			$qual=$2;
			$site=$1;
			next if ($qual<$minq);
			$nC++;
			if (!$seenC{$site}) {
				$seenC{$site}=1;
				$sC++;
			}
		}
		elsif ($_=~/\tchr14\t(\d+)\t(\d+)/) {
			$qual=$2;
			$site=$1;
			next if ($qual<$minq);
			$nD++;
			if (!$seenD{$site}) {
				$seenD{$site}=1;
				$sD++;
			}
		}
	}
	my $sumf=$nA+$nB+$nC+$nD;
	if ($sumf==0) { 
		print "$file\t$sH\t$sA\t$sB\t$sC\t$sD\t$nH\t$nA\t$nB\t$nC\t$nD\tNA\tNA\tNA\tNA\tNA\n";
		next;
	}
	my $fZ=sprintf("%.5f",$sumf/($sumf+$nH));
	$fA=sprintf("%.3f",$nA/$sumf);
	$fB=sprintf("%.3f",$nB/$sumf);
	$fC=sprintf("%.3f",$nC/$sumf);
	$fD=sprintf("%.3f",$nD/$sumf);
 	print "$file\t$sH\t$sA\t$sB\t$sC\t$sD\t$nH\t$nA\t$nB\t$nC\t$nD\t$fZ\t$fA\t$fB\t$fC\t$fD\n";
 }
 