#!/usr/bin/perl
#This program merges reads with the same primer and adds sample name to header of each read
use warnings;
use strict;


my @files = @ARGV;
my @filename_line;
my $path = "Output/Primer_trimming_sorting/" ;

foreach my $filename (@files){
    open FH, "$filename" or die;
    my $prefix = $filename ;
    $prefix =~ s/$path// ;
    $prefix =~ s/\.sorted.fasta//; 
#	@filename_line = split(/\./, $filename) ;
#	my $sample = $filename_line[0] ;
#	print "$sample\n" ;
    while (my $line = <FH>){
    chomp $line;
    	if ($line =~ s/@/>/) {
        	print "$line" . '-' . "$prefix\n";
   	 	} 
   	 	else {
        	print "$line\n";}
    	}
    	close FH;
}