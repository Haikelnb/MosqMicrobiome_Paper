#!/usr/bin/perl
#This program parses through the Final_output and only keeps unique species for the each of the queries best blast hit
use warnings;
use strict;


my $finaloutputfilename = $ARGV[0];
open (FINALOUTPUTFILE, $finaloutputfilename) or die "Can't open file '$finaloutputfilename' [$!]\n";


my %seen;
while (my $line = <FINALOUTPUTFILE>) {
    my $input_path = "Output/Final_Output/" ;
	my $output_path = "Output/Excel_Ready_Output/" ;
    my $prefix = $finaloutputfilename ;
    $prefix =~ s/$input_path//;
    $prefix =~ s/_Final.output.txt//;
	my $Excel_outfile = $output_path.$prefix."_Excel.Output.txt";
	open (FINALOUTFILE, '>>', $Excel_outfile) or die "Can't write to file '$Excel_outfile' [$!]\n";
    $seen{$line}++ or print FINALOUTFILE "$line";    ## This is the only line of code that matters in this script (smh) ##
}