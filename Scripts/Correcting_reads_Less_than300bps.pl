#!/usr/bin/perl
#This program makes sure paired reads R1 and R2 have the same amount of reads
use warnings;
use strict;

my $R1filename = $ARGV[0]; #assign a variable to file
my $R2filename = $ARGV[1];
#my $R1trimmedfilename = $ARGV[2];
#my $R2trimmedfilename = $ARGV[3];
open (R1FILE, $R1filename) ; #open file with filehandle
open (R2FILE, $R2filename);
#open (R1TRIMMEDFILE, '>', $R1trimmedfilename) or die "Can't write to file '$R1trimmedfilename' [$!]\n";    #opens filehandle for output file you-
#open (R2TRIMMEDFILE, '>', $R2trimmedfilename) or die "Can't write to file '$R2trimmedfilename' [$!]\n";    # -need to print stuff into


my @R1reads = <R1FILE>;
my @R2reads = <R2FILE>;


my @header_line2; 
my %header_hash2;
my $key2;

foreach (@R2reads) {
	   #every line of the Reads file is split by space and assigns a new variable for each line; the first line has two characters split by one space
	if (($_ =~ /HKL7MBCX2/) | ($_ =~ /M00708/)) {    #if the variable has the character 'HWI' 
		@header_line2 = split(" ", $_) ;
		$key2 = $header_line2[0];    #the first variable, which is the header of the reads and differentiates each read, is saved as a key
	} else {
		$header_hash2{$key2} .= $_ ; # the rest of the lines (2-4) are saved as values for the hash
	}
}
#foreach my $key (sort keys %header_hash2) { 
#	print "$key\n$header_hash2{$key}";          #can print the keys and values of each key
#}



my @header_line1;
my %header_hash1;
my $key1;

foreach (@R1reads) {
	if (($_ =~ /HKL7MBCX2/) | ($_ =~ /M00708/)) {
		@header_line1 = split(" ", $_); #every line of the Reads file is split by space; this splits them into two characters
		$key1 = $header_line1[0];		# the first element which is the part of the header of each sequence that is unique is saved into a varialbe
	} else {
		$header_hash1{$key1} .= $_ ;	#saves every line until the if statment is true again into the hash with the given key
	}
}



foreach my $key (sort keys %header_hash1) {
	my $prefix = $R1filename ; 
   	$prefix =~ s/\.fastq//;
	my $R1trimmedfilename = $prefix.".corrected.fastq";		# this is for naming the outfile
	if (exists $header_hash2{$key}  ) {   #if the key for R2 reads exists in the keys of R1, the R1s are printed in the file
		open (R1TRIMMEDFILE, '>>', $R1trimmedfilename) or die "Can't write to file '$R1trimmedfilename' [$!]\n";
		print R1TRIMMEDFILE "$key\n$header_hash1{$key}";
			
	}	
}

foreach my $key (sort keys %header_hash2) {
	my $prefix = $R2filename ; 
   	$prefix =~ s/\.fastq//;
	my $R2trimmedfilename = $prefix.".corrected.fastq";
	if (exists $header_hash1{$key} ) {     #if the key for R1 reads exists in the keys of R2, the R2s are printed in the file
		open (R2TRIMMEDFILE, '>>', $R2trimmedfilename) or die "Can't write to file '$R2trimmedfilename' [$!]\n"; 
		print R2TRIMMEDFILE "$key\n$header_hash2{$key}";   
			
	} 
}

close R1TRIMMEDFILE;
close R2TRIMMEDFILE;
close R1FILE;
close R2FILE;

exit;