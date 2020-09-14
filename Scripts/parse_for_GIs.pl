#!/usr/bin/perl
#This program parses through the blast_output file (just like the parse_blastOutput script) and prints out the gi_ids only in a single column
use warnings;					## The comments for this script can be find in parse_blastOutput script with minor changes ## 
use strict;


my $blastfilename = $ARGV[0];
open (BLASTFILE, $blastfilename) or die "Can't open file '$blastfilename' [$!]\n";

my $query;
my $highest_rawscore = 0;

while (my $line = <BLASTFILE>) {
	 my $blast_outfile = $blastfilename.".parsed";
	 chomp $line;
	 if ($line =~ /Query/) {
	 	$query = substr ($line, 9) ;	 	
	 }  
	 if (defined ($query) ) { 
		unless ($line =~ /#/ ) {
	 		my @blast_line = split("\t", $line);
	 		my $rawscore = $blast_line[6] ;
	 		my $gi_id = $blast_line[1] ;
	 		
			if ($highest_rawscore == 0 ) {
				$highest_rawscore = $rawscore ;
			}
	 		
	 		my @GI_trimmed;					# This array is initialized for extracting only the GI_IDs
	 		if ($rawscore >= $highest_rawscore) {
	 			@GI_trimmed = split(/\|/, $gi_id) ;
	 			my $input_path = "Output/Blast_outputs/" ;
				my $output_path = "Output/GI_IDs/" ;  
	 			my $prefix = $blastfilename ; 
    			$prefix =~ s/\.blast_output//;
    			$prefix =~ s/$input_path//;
	 			my $giIDs_outfile = $output_path.$prefix."_gi_id.txt";		# this is for naming the outfile
	 			open (OUTFILE, '>>', $giIDs_outfile) or die "Can't write to file '$giIDs_outfile' [$!]\n";
	 			print OUTFILE "$GI_trimmed[1]\n" ;
	
				$highest_rawscore = $rawscore ;
			}		 
		}

	 } 
	 if ($line =~ /#/ ) {
	 	$highest_rawscore = 0 ;
	 }	 	 
}

close OUTFILE ;
close BLASTFILE ;