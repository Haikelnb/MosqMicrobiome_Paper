#!/usr/bin/perl
#This program parses through each line of Results file and determines if the sequences are mutant/WT sequences and prints out the genotype per sample
use warnings;
use strict;


my $resultsfilename = $ARGV[0]; 
open (RESULTSFILE, $resultsfilename) or die "Can't open file '$resultsfilename' [$!]\n";
my @sample = <RESULTSFILE> ;



my $KDR_WT = "ACGACTAAATTTC" ;
my $kdr_w = "ACGACAAAATTTC" ;
my $kdr_e = "ACGACTGAATTTC" ;

my $prefix = $resultsfilename ;
$prefix =~ s/\.txt// ;

foreach (@sample) {
	my @sample_line = split("\t", $_) ;
	if ($sample_line[3] =~ /$KDR_WT/) {
		print "\nKDR_WT" ;
	}
	elsif ($sample_line[3] =~ /$kdr_w/) {
		print "\nkdr_w" ;
	}
	elsif ($sample_line[3] =~ /$kdr_e/) {
	  	print "\nkdr_e" ; 
	}
	else {
	  	print "\nother" ;
	}
	if ( defined($sample_line[5])) {
		if ($sample_line[5] =~ /$KDR_WT/) {
			print "\tKDR_WT\n" ;
		}
		elsif ($sample_line[5] =~ /$kdr_w/) {
			print "\tkdr_w\n" ;
		}
		elsif ($sample_line[5] =~ /$kdr_e/) {
	  		print "\tkdr_e\n" ; 
		}
		else {
	  		print "\tother\n" ;
		}
	}
}
