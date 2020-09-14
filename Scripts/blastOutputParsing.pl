#!/usr/bin/perl
#This program parses through the unique_sequences blast output and counts # of samples for each query sequence
use warnings;
use strict;


my $blastfilename = $ARGV[0];
my $duplicatesfilname = $ARGV[1]; 
open (BLASTFILE, $blastfilename) or die "Can't open file '$blastfilename' [$!]\n";
open (DUPLICATESFILE, $duplicatesfilname) or die "Can't open file '$duplicatesfilname' [$!]\n";


my %header_hash;
my $header_name;
my $sample_name;
my $header_sample_name; 
while (my $duplicate_line = <DUPLICATESFILE>) {
	chomp $duplicate_line;
	my @duplicate_headers = split(/\./, $duplicate_line) ;  ## need to be careful with the "1-" splitter, all reads files should have this to separate the header and sample name
#	print "$duplicate_headers[1]\n" ;
	if ( $duplicate_line =~ /^[>]/ ) {
	  	$header_name = $duplicate_line ;
	  	$header_name =~ s/>// ;	
	  	$header_sample_name = $duplicate_headers[1] ;
	  	$sample_name = $header_sample_name ;
	  	$header_hash{$header_name}{$sample_name}++ ;
	}
	if ( $duplicate_line =~ /^[+]/ ) {
		$sample_name = $duplicate_headers[1] ;
	  	$header_hash{$header_name}{$sample_name}++ ;	
	}
}

#foreach my $header_name (keys %header_hash) {
#	foreach my $name (keys %{$header_hash{$header_name}}) {	
#		print "$name\t$header_name\tcount:$header_hash{$header_name}{$name}\n" ;
#	}
#}
		
		
my $query;
my $highest_rawscore = 0;

while (my $line = <BLASTFILE>) {
#	 my $blast_outfile = $blastfilename.".parsed";
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
	 		
			if ($rawscore >= $highest_rawscore) {
				foreach my $name (keys %{ $header_hash{$query} }) {
					my $prefix = $blastfilename ;
					my $blast_outfile = $prefix.".parsed";   # the new name for the outfile
					open (BLASTOUTFILE, '>>', $blast_outfile) or die "Can't write to file '$blast_outfile' [$!]\n";	
					print BLASTOUTFILE "$name\tcount:$header_hash{$query}{$name}\t\t$line\n" ;
#					print "$name\tcount:$header_hash{$query}{$name}\t\t$line\n" ;
				}
				$highest_rawscore = $rawscore ;
			}
		}
	 } 
	if ($line =~ /#/ ) {
	 	$highest_rawscore = 0 ;
	}	 	 
}

close BLASTOUTFILE;
close BLASTFILE;
close DUPLICATESFILE; 