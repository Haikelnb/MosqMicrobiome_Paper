#!/usr/bin/perl
#This program counts the to identity the top (highest count) 2 sequences in each sample 
use warnings;										
use strict;

my $sortedfilename = $ARGV[0]; 
open (SORTEDFILE, $sortedfilename) or die "Can't open file '$sortedfilename' [$!]\n";

my %sequence_hash;
my $header;
my $sequence_line;

my $Total_Reads_count; 
my %Reads_count;

while (my $file_line = <SORTEDFILE>) {
	chomp $file_line;
	unless ($file_line =~ /^@/) {
		$Total_Reads_count++ ;    # counts the number of reads
		$sequence_line = $file_line ;        # the lines that are not headers are saved as sequence
#		$sequence_hash{$sequence_line} = $header;    # this is how you normally save a hash
#		push ( @{ $sequence_hash{$sequence_line} }, $sequence_line );    # this line creates array of hashes for the same sequences that have multiple
		$sequence_hash{$sequence_line}++ ;												# headers. Also the key is the sequence and the headers are values
	}															  
}

my $line_count = 1;
my $seq1;
my $seq2;
my $seq1_sequence;
my $seq2_sequence;
foreach my $key (sort  { $sequence_hash{$b} <=> $sequence_hash{$a} } keys %sequence_hash) {
	if ($Total_Reads_count > 100) {
		if ($line_count == 1) {
#			print "$line_count\t$key\t$sequence_hash{$key}\t$Total_Reads_count\n";
			$seq1 = $sequence_hash{$key} ;
			$seq1_sequence = $key ;
			$line_count++ ;
			
		} 
		elsif ($line_count == 2) {
#			print "$line_count\t$key\t$sequence_hash{$key}\t$Total_Reads_count\n";
			$seq2 = $sequence_hash{$key} ;
			$seq2_sequence = $key ;
			$line_count++ ;
		} 
	} 
}

my $prefix = $sortedfilename ;
$prefix =~ s/\.sorted.fasta// ;
if (defined ($seq1 and $seq2)) {
	my $cutoff_calcuation = ($seq2/($seq1+$seq2)) ;
#	print "$prefix\t$cutoff_calcuation\n" ;
	if ($cutoff_calcuation <= 0.25) {
		print "$prefix\tTotal:$Total_Reads_count\tSeq1:$seq1\t$seq1_sequence\n" ;
	}
	if ($cutoff_calcuation > 0.25) {
		print "$prefix\tTotal:$Total_Reads_count\tSeq1:$seq1\t$seq1_sequence\tSeq2:$seq2\t$seq2_sequence\n" ;
	}
}