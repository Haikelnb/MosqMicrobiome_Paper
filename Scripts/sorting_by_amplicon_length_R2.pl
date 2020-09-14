#!/usr/bin/perl
# This program counts the # of reads per primer (R1s) from the raw_reads file before demultiplexing and marks where the second primer is found
use warnings;
use strict;

my $primerfilename = $ARGV[0];
my $rawreadsfilename = $ARGV[1];
open (RAWREADSFILE, $rawreadsfilename) ;     
open (PRIMERFILE, $primerfilename) ;	 # opens the file with all the primers (names and sequences)
my @reads = <RAWREADSFILE> ;
my @primers = <PRIMERFILE> ;

my %Primer_hash;
my @primer_line;
my $primer_key;
my $primerR;
my $primerF;
my $amplicon_length;

foreach (@primers) {
	@primer_line = split("\t", $_); ####### have to make sure that primers in txt file are saperated by tabs, there should a tab after the last primer for each line ######

	$primerF = substr ($primer_line[2], 0, 10) ;
	$primerR = substr ($primer_line[3], 0, 10) ;
	$primerR = IUB_to_regexp($primerR) ;
	$primerF = IUB_to_regexp($primerF) ;
	
		
	$primer_key = $primer_line[0];     # the first element (of each line) which is the primer name is saved as primer_key
	$amplicon_length = $primer_line[1] ;
	
#	$Primer_hash{$primer_key}{Reverse_primer} = $primerR ;  		# the reverse primer (the sequence needed for matching) is saved 
#	$Primer_hash{$primer_key}{Amplicon_length} = $amplicon_length;
	$Primer_hash{$primer_key}{Forward_primer} = $primerF; 
	$Primer_hash{$primer_key}{Amplicon_length} = $amplicon_length;    
}

#foreach my $primer_key (keys %Primer_hash) {
#	foreach my $amplicon_length (keys %{$Primer_hash{$primer_key}}) {	
#		print "$primer_key\t$Primer_hash{$primer_key}{Amplicon_length}\n" ;
#	}
#}


#foreach my $primer_key (keys %Primer_hash) {
#	foreach my $amplicon_length (keys %{$Primer_hash{$primer_key}}) {	
#		if ( $amplicon_length = 350 ) {
#			print "$primer_key\t$amplicon_length\n" ;
#		}
#	}
#}

my $trimmed_sequence;
my $trimmed_quality;
my $untrimmed_read;

my $trimmed_hit_count;
my $untrimmed_hit_count;


my $header_line = "";			# this line is needed to capture the line above
my $count = 0;	
my $count_2 = 0;
for (my $i = 0; $i < @reads; $i++) {				# loop through file line by line
	my $prefix = $rawreadsfilename ; 
#   	my $path = $ARGV[2];
   	$prefix =~ s/\_Untrimmed.fastq//;
	my $trimmed_R2_file = $prefix."_trimmed_50bps.fastq";		# this is for naming the outfile
	my $Untrimmed_R2_file = $prefix."_greaterThan300_ampl.fastq";
	foreach my $primer_key (sort keys %Primer_hash) {
#		print "$Primer_hash{$key}{Reverse_primer}\t$Primer_hash{$key}{Forward_primer}\n";
	 	if ( ($Primer_hash{$primer_key}{Amplicon_length} =~ /Less_than_300/)  && $reads[$i] =~ /^$Primer_hash{$primer_key}{Forward_primer}/ ) {    # if the read matches the reverse primer and the amplicon is shorter than 300bps	
			$count = 0 ;
			
			$trimmed_sequence = substr($reads[$i], 0, -51) ; 
			open (TRIMMEDOUTFILE, '>>', $trimmed_R2_file) or die "Can't write to file '$trimmed_R2_file' [$!]\n";
			print TRIMMEDOUTFILE "$header_line";
			print TRIMMEDOUTFILE "$trimmed_sequence\n" ;
			
			$count++ ;
			++$trimmed_hit_count ;
		}
		if ( ($Primer_hash{$primer_key}{Amplicon_length} =~ /Greater_than_300/)  && $reads[$i] =~ /^$Primer_hash{$primer_key}{Forward_primer}/ ) {    # if the read matches the reverse primer and the amplicon is longer than 300bps	
			$count_2 = 0 ;
			
			open (UNTRIMMEDOUTFILE, '>>', $Untrimmed_R2_file) or die "Can't write to file '$Untrimmed_R2_file' [$!]\n";
			print UNTRIMMEDOUTFILE "$header_line" ;
			$untrimmed_read = $reads[$i] ;
			print UNTRIMMEDOUTFILE "$untrimmed_read" ;
			
			$count_2++ ;
			++$untrimmed_hit_count ;
		}
	}
	
	if ( defined($trimmed_sequence)  && $count == 2) {    
		print TRIMMEDOUTFILE "$reads[$i]" ;
	}
	elsif ( defined($untrimmed_read) && $count_2 == 2) {  
		print UNTRIMMEDOUTFILE "$reads[$i]" ;
	}
	if ( defined($trimmed_sequence)  && $count == 3) {
		$trimmed_quality = substr($reads[$i], 0, -51) ;   # this is for trimming the quality line the same as the read
		print TRIMMEDOUTFILE "$trimmed_quality\n" ;
	}
	elsif ( defined($untrimmed_read) && $count_2 == 3) {  
		print UNTRIMMEDOUTFILE "$reads[$i]" ;
	}
	
	$header_line = $reads[$i];
	$count++ ; 
	$count_2++ ;
}

if ( defined ($trimmed_hit_count) and defined ($untrimmed_hit_count)) {
	print "$rawreadsfilename\t50bps Trimmed Reads Hit Count:$trimmed_hit_count\tUntrimmed Reads Hit Count:$untrimmed_hit_count\n" ;
}
elsif ( !defined ($trimmed_hit_count) ) {
	print "$rawreadsfilename\t50bps Trimmed Reads Hit Count:0\tUntrimmed Reads Hit Count:$untrimmed_hit_count\n" ;
}
elsif ( !defined ($untrimmed_hit_count) ) {
	print "$rawreadsfilename\t50bps Trimmed Reads Hit Count:$trimmed_hit_count\tUntrimmed Reads Hit Count:0\n" ;
}

sub IUB_to_regexp {
	my($iub) = @_;
	my $regular_expression = '';
	my %iub2character_class = (
		A => 'A',
		C => 'C',
		G => 'G',
		T => 'T',
		R => '[GA]',
		Y => '[CT]',
		M => '[AC]',
		S => '[GC]',
		W => '[AT]',
		B => '[CGT]',
		D => '[AGT]',
		H => '[ACT]',
		K => '[GT]',
		V => '[ACG]',
		N => '[ACGT]',
	);
	for (my $i = 0; $i < length($iub); ++$i) {
		$regular_expression .= $iub2character_class{substr($iub, $i, 1)} ;
	}
	return $regular_expression ;
}
	
close TRIMMEDOUTFILE; 
close UNTRIMMEDOUTFILE;
close RAWREADSFILE ;
close PRIMERFILE ;