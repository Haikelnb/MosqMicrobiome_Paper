#!/usr/bin/perl
# This program counts the # of reads per primer from the raw_reads file before demultiplexing
use warnings;
use strict;

my $primerfilename = $ARGV[0];
my $mergedfilename = $ARGV[1];
open (MERGEDFILE, $mergedfilename) ;     # opens the merged sequences file
open (PRIMERFILE, $primerfilename) ;	 # opens the file with all the primers (names and sequences
my @reads = <MERGEDFILE> ;
my @primers = <PRIMERFILE> ;


my %Primer_hash;
my @primer_line;
my $primer_key;
my $primerR;
my $primerF;

foreach (@primers) {
	@primer_line = split("\t", $_); ####### have to make sure that primers in txt file are saperated by tabs, there should a tab after the last primer for each line ######

	$primerR = substr ($primer_line[1], 0, 10) ;
	$primerF = substr ($primer_line[2], 0, 10) ;
	$primerR = IUB_to_regexp($primerR) ;
	$primerF = IUB_to_regexp($primerF) ;
		
	$primer_key = $primer_line[0];     # the first element (of each line) which is the primer name is saved as primer_key
	$Primer_hash{$primer_key}{Reverse_primer} = $primerR; 		# the reverse primer (the sequence needed for matching) is saved 
	$Primer_hash{$primer_key}{Reverse_primer_for_length} = $primer_line[1];
	$Primer_hash{$primer_key}{Forward_primer} = $primerF; 
	$Primer_hash{$primer_key}{Forward_primer_for_length} = $primer_line[2];    
}

#foreach my $key (sort keys %Primer_hash) {
#	print "$key\t$Primer_hash{$key}{Reverse_primer_for_length}\t$Primer_hash{$key}{Forward_primer_for_length}\n"
#}

my %Primer_count_hash;

my $trimmed_sequence;
my $Forward_primer_length;
my $Reverse_primer_length;

my $hit_outfile;
my $path;
my $without_path;
my $header_line = "";								# this line is needed to capture the line above
for (my $i = 0; $i < @reads; $i++) {				# loop through file line by line
	foreach my $key (sort keys %Primer_hash) {
#	print "$Primer_hash{$key}{Reverse_primer}\t$Primer_hash{$key}{Forward_primer}\n";
	 	if ($reads[$i] =~ /^$Primer_hash{$key}{Reverse_primer}/  && $reads[$i] =~ /$Primer_hash{$key}{Forward_primer}/ ) {    # both the reverse and forward primers must match
	 		$Reverse_primer_length = length $Primer_hash{$key}{Reverse_primer_for_length};   # this and the next line is for trimming the primer sequences from the matched sequence
	 		$Forward_primer_length = length $Primer_hash{$key}{Forward_primer_for_length};
#			$trimmed_sequence = $reads[$i] =~ /$Primer_hash{$key}{Reverse_primer}(.*)$Primer_hash{$key}{Forward_primer}/ ;   this saves the sequence between two regexp matches
			$trimmed_sequence = substr($reads[$i], $Reverse_primer_length, -($Forward_primer_length + 1)) ;  # this captures the sequence in between the primers
	 		
			my $sequence_length = length ($trimmed_sequence) ;
			if ($sequence_length > 90) {
				my $prefix = $mergedfilename ; 
	 			$path = "Output/Primer_trimming_sorting/" ;   #Output/Primer_trimming_sorting/
	 			$prefix =~ s/$path// ;
   				$prefix =~ s/\.extendedFrags.fastq// ;
   				$without_path = $key.".$prefix.fasta" ;
	 			$hit_outfile = $path.$key.".$prefix.sorted.fasta" ;		# this is for naming the outfile
	 			open (OUTFILE, '>>', $hit_outfile) or die "Can't write to file '$hit_outfile' [$!]\n";
				print OUTFILE "$header_line";    # prints the header followed by the trimmed matched sequence
				print OUTFILE "$trimmed_sequence\n";
				
				$Primer_count_hash{$without_path}{Real_hit}++ ;

			}
			elsif ($sequence_length <= 90)  {
				my $prefix = $mergedfilename ; 
	 			$path = "Output/Primer_trimming_sorting/" ;   #Output/Primer_trimming_sorting/
	 			$prefix =~ s/$path// ;
   				$prefix =~ s/\.extendedFrags.fastq// ;
   				$without_path = $key.".$prefix.fasta" ;
				$Primer_count_hash{$without_path}{Short_hit}++ ;
			}	
		}		
	}
	$header_line = $reads[$i];  # saves the current line as the header line ( therefore when we match the sequence line the line above is always the header)
#	++$total_lines_count;		# counts every line in the file
}



#my $Reads_NoHits = (($total_lines_count/4) - $primer_hit_count) ;  #total num. of lines is divided by 4 b/c every read has four lines
#my $total_reads = ($total_lines_count/4) ;
#print "$mergedfilename Number of sequences with >=1 Hits $primer_hit_count\n";
#print "$mergedfilename Number of sequences with 0 Hits $Reads_NoHits\tTotal:$total_reads\n";

foreach my $without_path (sort keys %Primer_count_hash) {
	if ( defined ($Primer_count_hash{$without_path}{Real_hit}) and defined ($Primer_count_hash{$without_path}{Short_hit}) ) {
		print "$without_path\tReal_Hits: $Primer_count_hash{$without_path}{Real_hit}\tShort_hits: $Primer_count_hash{$without_path}{Short_hit}\n" ;
	}
	elsif ( defined ($Primer_count_hash{$without_path}{Real_hit}) and !defined($Primer_count_hash{$without_path}{Short_hit}) ) {
		print "$without_path\tReal_Hits: $Primer_count_hash{$without_path}{Real_hit}\tShort_hits: 0\n" ;
	}
	elsif ( !defined($Primer_count_hash{$without_path}{Real_hit}) and defined($Primer_count_hash{$without_path}{Short_hit}) ) {
		print "$without_path\tReal_Hits: 0\tShort_hits: $Primer_count_hash{$without_path}{Short_hit}\n" ;
	}
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


close OUTFILE ;
close MERGEDFILE ;
close PRIMERFILE ;
