#!/usr/bin/perl
#This program keeps only unique sequences from a fasta file and counts the number the sequences (also replaces @ with >), 
use warnings;										#(AND saves all the sequences with multiple headers)
use strict;

my $Fastafilename = $ARGV[0];
open (FASTAFILE, $Fastafilename) or die "Can't write to file '$Fastafilename' [$!]\n";
#my @sequences = <FASTAFILE> ;


my %sequence_hash;
my $header;
my $sequence_line;

my $Total_Reads_count; 

my %count;

while (my $line = <FASTAFILE>){
    chomp $line;
	if ($line =~ />/) {   #if the lines has an @ it replaces it with > (to make it a fasta file)
		$header = $line ;         # this line is saved as a header
		$Total_Reads_count++ ;    # counts the number of reads
	}
	else {
		$sequence_line = $line ;        # the lines that are not headers are saved as sequence
#		$sequence_hash{$sequence_line} = $header;    # this is how you normally save a hash
		push ( @{ $sequence_hash{$sequence_line} }, $header );    # this line creates array of hashes for the same sequences that have multiple
		$count{$sequence_hash{$sequence_line}}++ ;
	}															  # headers. Also the key is the sequence and the headers are values
}


my $counter = 0 ;           # for counting unique sequences

my $prefix = $Fastafilename ;  
foreach my $sequence_line (sort keys %sequence_hash) {
#	if ($count{$sequence_hash{$sequence_line}} >= 10 ) {
#		my $prefix = $Fastafilename ;     		# the prefix is needed to keep bases of input file name for the output file name
		$prefix =~ s/\.AHBGM_20180417_M00708_IL100101409_S1_L001.fasta//;       # removes the fasta from the file name so it can be renamed
		my $unique_sequence_outfile = $prefix.".unique.fasta";   # the new name for the outfile
		open (UNIQUEOUTFILE, '>>', $unique_sequence_outfile) or die "Can't write to file '$unique_sequence_outfile' [$!]\n";
#		print "@{ $sequence_hash{$sequence_line} }[0]\t$count{$sequence_hash{$sequence_line}}\n$sequence_line\n" ;
		print UNIQUEOUTFILE "@{ $sequence_hash{$sequence_line} }[0]\t$count{$sequence_hash{$sequence_line}}\n$sequence_line\n" ;
		$counter++ ; 
#	}		
}

print "$prefix\tUnique-sequences $counter\tTotal-sequences $Total_Reads_count\n" ;

close UNIQUEOUTFILE ;
close FASTAFILE;

