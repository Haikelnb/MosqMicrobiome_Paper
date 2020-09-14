#!/usr/bin/perl
# This program trims R2 reads that have the second (Reverse) primer when it is present; 
use warnings;        
use strict;

#### The primer_file for R2_reads needs to have primer name, forward primer, then Reverse Complement of Reverse Primer, tab delimited ####
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

foreach (@primers) {
	@primer_line = split("\t", $_); ####### have to make sure that primers in txt file are saperated by tabs, there should a tab after the last primer for each line ######

	$primerF = substr ($primer_line[1], 0, 10) ;
	$primerR = substr ($primer_line[2], 0, 10) ;
	$primerF = IUB_to_regexp($primerF) ;
	$primerR = IUB_to_regexp($primerR) ;
		
	$primer_key = $primer_line[0];     # the first element (of each line) which is the primer name is saved as primer_key
	$Primer_hash{$primer_key}{Forward_primer} = $primerF; 		# the reverse primer (the sequence needed for matching) is saved 
#	$Primer_hash{$primer_key}{Forward_primer_for_length} = $primer_line[1];
	$Primer_hash{$primer_key}{Reverse_primer} = $primerR; 
#	$Primer_hash{$primer_key}{Reverse_primer_for_length} = $primer_line[2];    
}


my $trimmed_sequence;
my $trimmed_quality;
my $untrimmed_read;
my $position;

my $trimmed_hit_count;
my $untrimmed_hit_count; 
my $First_primer_hit;

my $header_line = "";			# this line is needed to capture the line above
my $count = 0;
my $count_2 = 0;								
for (my $i = 0; $i < @reads; $i++) {				# loop through file line by line
	my $prefix = $rawreadsfilename ; 
	my $path = $ARGV[2];
   	$prefix =~ s/\.fastq//;
	my $trimmed_R2_file = $prefix."_trimmed_secondPrimer.fastq";		# this is for naming the outfile
	my $Untrimmed_R2_file = $prefix."_Untrimmed.fastq";
	foreach my $key (sort keys %Primer_hash) {
#	print "$Primer_hash{$key}{Reverse_primer}\t$Primer_hash{$key}{Forward_primer}\n";
	 	if ($reads[$i] =~ /^$Primer_hash{$key}{Forward_primer}/  && $reads[$i] =~ /$Primer_hash{$key}{Reverse_primer}/ ) {    # both the reverse and forward primers must match
			$count = 0 ;
#			$trimmed_sequence = $reads[$i] =~ /$Primer_hash{$key}{Reverse_primer}(.*)$Primer_hash{$key}{Forward_primer}/ ;   this saves the sequence between two regexp matches
			$position = $+[0] ;   #marks the position of the last successful match 
#			print "$position\n" ;
			$trimmed_sequence = substr($reads[$i], 0, $position) ;  # this captures the sequence before the second primer
			open (TRIMMEDOUTFILE, '>>', $path . "Trimmed/" . $trimmed_R2_file) or die "Can't write to file '$trimmed_R2_file' [$!]\n";
			print TRIMMEDOUTFILE "$header_line";
#			print "$reads[$i]" ;
			print TRIMMEDOUTFILE "$trimmed_sequence\n" ;
			
			$count++ ;
			++$trimmed_hit_count ;
		}
		if 	($reads[$i] =~ /^$Primer_hash{$key}{Forward_primer}/ ) {
			++$First_primer_hit;
			unless ($reads[$i] =~ /$Primer_hash{$key}{Reverse_primer}/ ) { 	
				$count_2 = 0 ;
				open (UNTRIMMEDOUTFILE, '>>', $path . "Untrimmed/" . $Untrimmed_R2_file) or die "Can't write to file '$Untrimmed_R2_file' [$!]\n";
				print UNTRIMMEDOUTFILE "$header_line";
				$untrimmed_read = $reads[$i] ;
				print UNTRIMMEDOUTFILE "$untrimmed_read" ;
				$count_2++ ;
				++$untrimmed_hit_count ;
			}
		}	
	}
	if ( defined($trimmed_sequence)  && $count == 2) {    
		print TRIMMEDOUTFILE "$reads[$i]" ;
	}
	elsif ( defined($untrimmed_read) && $count_2 == 2) {  
		print UNTRIMMEDOUTFILE "$reads[$i]" ;
	}
	if ( defined($trimmed_sequence)  && $count == 3) {
		$trimmed_quality = substr($reads[$i], 0, $position) ;   # this is for trimming the quality line the same as the read
		print TRIMMEDOUTFILE "$trimmed_quality\n" ;
	}
	elsif ( defined($untrimmed_read) && $count_2 == 3) {  
		print UNTRIMMEDOUTFILE "$reads[$i]" ;
	}
	
	$header_line = $reads[$i];  # saves the current line as the header line ( therefore when we match the sequence line the line above is always the header)
	$count++ ;                   # this counts line so the second and thrid line of a read (+ sign and quality) is printed for each read
	$count_2++ ;
}

if ( defined($First_primer_hit and $trimmed_hit_count and $untrimmed_hit_count ) ) {
	print "$rawreadsfilename\tTotal_First_Primer_Hit:$First_primer_hit\tTrimmed Reads Hit Count:$trimmed_hit_count\tUntrimmed Reads Hit Count:$untrimmed_hit_count\n" ;
}
elsif ( !defined ($First_primer_hit)) {
	print "$rawreadsfilename\tTotal_First_Primer_Hit:0\tTrimmed Reads Hit Count:0\tUntrimmed Reads Hit Count:0\n" ;
}
elsif ( !defined ($trimmed_hit_count)) {
	print "$rawreadsfilename\tTotal_First_Primer_Hit:$First_primer_hit\tTrimmed Reads Hit Count:0\tUntrimmed Reads Hit Count:$untrimmed_hit_count\n" ;
}
elsif ( !defined ($untrimmed_hit_count)) {
	print "$rawreadsfilename\tTotal_First_Primer_Hit:$First_primer_hit\tTrimmed Reads Hit Count:$trimmed_hit_count\tUntrimmed Reads Hit Count:0\n" ;
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

