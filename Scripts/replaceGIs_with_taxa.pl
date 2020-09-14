#!/usr/bin/perl
#This program parses through the unique_sequences blast output and counts # of samples for each query sequence
use warnings;
use strict;


my $parsedfilename = $ARGV[0];
my $taxafilename = $ARGV[1]; 
open (PARSEDFILE, $parsedfilename) or die "Can't open file '$parsedfilename' [$!]\n";
open (TAXAFILE, $taxafilename) or die "Can't open file '$taxafilename' [$!]\n";


my %giID_hash;
my $gi_id ;

while (my $line = <TAXAFILE>) {
	chomp $line ;
	unless ($line =~ /^gi/) {
#		print "$line\n" ;
		my @taxa_line = split("\t", $line);
		unless ($taxa_line[1] =~ /sp./ or $taxa_line[1] =~ /cf./ or $taxa_line[1] =~ /str./ ) { 
			$gi_id = $taxa_line[0] ;
			my $taxa_data = join "\t", @taxa_line[1..9];
			$giID_hash{$gi_id} = $taxa_data ;
		}
	}
}




while (my $blast_line = <PARSEDFILE>) {
	chomp $blast_line ;
	my @parsed_line = split(" ", $blast_line) ;
	foreach $gi_id ( keys %giID_hash) {
		if ($blast_line =~ /$gi_id/ ) {
			my $prefix = $parsedfilename ;
    		$prefix =~ s/\.blast_output.parsed//;
	 		my $Final_outfile = $prefix."_Final.output.txt";		# this is for naming the outfile
	 		open (FINALOUTFILE, '>>', $Final_outfile) or die "Can't write to file '$Final_outfile' [$!]\n";
			print FINALOUTFILE "$parsed_line[0] $parsed_line[1]\t$parsed_line[2]\t$giID_hash{$gi_id}\t$parsed_line[4]\t$parsed_line[5]\t$parsed_line[6]\t$parsed_line[7]\t$parsed_line[8]\n" ;
		
		}
	}
}

close FINALOUTFILE ;
close PARSEDFILE ;
close TAXAFILE ;


 


	