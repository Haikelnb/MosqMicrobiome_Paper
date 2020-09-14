#!/usr/bin/perl
#This program keeps only unique sequences from a fasta file and counts the number the sequences (also replaces @ with >), 
use warnings;										#(AND saves all the sequences with multiple headers)
use strict;

my $Fastafilename = $ARGV[0];
open (FASTAFILE, $Fastafilename) or die "Can't write to file '$Fastafilename' [$!]\n";
my @sequences = <FASTAFILE> ;


my %sequence_hash;
my $header;
my $sequence_line;
my $input_path;
my $output_path;
my $without_path;


#foreach(my $line = <FASTAFILE>) {
#	chomp $line;
my $Total_Reads_count; 
my %Reads_count;

foreach (@sequences) {
	chomp;
	if ($_ =~ />/) {   #if the lines has an @ it replaces it with > (to make it a fasta file)
		$header = $_ ;         # this line is saved as a header
		$Total_Reads_count++ ;    # counts the number of reads
	}
	else {
		$sequence_line = $_ ;        # the lines that are not headers are saved as sequence
#		$sequence_hash{$sequence_line} = $header;    # this is how you normally save a hash
		push ( @{ $sequence_hash{$sequence_line} }, $header );    # this line creates array of hashes for the same sequences that have multiple
		$Reads_count{$sequence_hash{$sequence_line}}++ ;												# headers. Also the key is the sequence and the headers are values
	}															  
}

 
my $count = 0 ;         # for counting unique sequences


foreach my $sequence_line (sort keys %sequence_hash) {		# goes through each key which are the sequences
	if ($Reads_count{$sequence_hash{$sequence_line}} >= 10 ) {
		$input_path = "Output/Merged_by_Primer/" ;
		$output_path = "Output/Unique_Sequences/" ;
		my $prefix = $Fastafilename ;     		# the prefix is needed to keep bases of input file name for the output file name
		$prefix =~ s/$input_path// ;
		$prefix =~ s/\.sorted.merged.fasta//;       # removes the fasta from the file name so it can be renamed
		$without_path = $prefix.".unique.fasta";
		my $unique_sequence_outfile = $output_path.$prefix.".unique.fasta";   # the new name for the outfile
		open (UNIQUEOUTFILE, '>>', $unique_sequence_outfile) or die "Can't write to file '$unique_sequence_outfile' [$!]\n";
		print UNIQUEOUTFILE "@{ $sequence_hash{$sequence_line} }[0]\t$Reads_count{$sequence_hash{$sequence_line}}\n" ;     # the [0] prints the first value only (first header)
		print UNIQUEOUTFILE "$sequence_line\n" ;
#		print "@{ $sequence_hash{$sequence_line} }[0]$sequence_line" ; "@{ $sequence_hash{$sequence_line} }[0]\t$Reads_count{$sequence_hash{$sequence_line}}\n$sequence_line\n" ;
	
	
		my $mulitple_header_outfile = $output_path.$prefix.".duplicate.fasta";
		open (DUPLICATE, '>>', $mulitple_header_outfile) or die "Can't write to file '$mulitple_header_outfile' [$!]\n";
#		print DUPLICATE ( @{ $sequence_hash{$sequence_line} }, "$sequence_line"  ) ;
 		print DUPLICATE  join ( "\n+", @{ $sequence_hash{$sequence_line} } ) ;    # this prints all the headers and the ones after the first
 		print DUPLICATE "\n$sequence_line\n"   ;	#prints the sequence after the header(s)	# header start with a "+" to find them quick in afile
 		$count++ ;            #counts occurrence of each unique sequence (every repeated sequence is represented once b/c they are the keys of the hash
	}
}

if ($count > 0) {
	print "$without_path\tUnique-sequences $count\tTotal-sequences $Total_Reads_count\n" ;
}
else {
	print "$Fastafilename\tUnique-sequences 0\tTotal-sequences $Total_Reads_count\n" ;
}

close UNIQUEOUTFILE;
close DUPLICATE;
close FASTAFILE;
