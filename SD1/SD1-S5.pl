#convert a gene list GFF3 file to the trimmed format used by COMPILE

#for these studies, converts maize transcript names to gene names on line 42 by removing _T suffix
#feature IDs to retain specified on lines 32 and 38
#input file name specified on line 16

#STRUCTURE
#SD1/
#├── SD1-S5.pl #script to run to perform GFF cleaning
#├── input .gff3 file
#└── GeneNameList.gff3 output file

use warnings;
use strict;
	
open my $input, '<', 'Zea_mays.B73_RefGen_v4.45.gff3'; #change filename as necessary
open my $output, '>', "GeneNameList.gff3";

while (<$input>) {

	next unless($_ !~ /^#/);

	chomp $_;
	my @t = split("\t", $_);
	my $chromosome = $t[0];
	my $source = $t[1];
	my $type = $t[2];
	my $start = $t[3];
	my $stop = $t[4];
	my $id = $t[8];
	
	if (($type =~ m/snRNA|snoRNA|lnc_RNA|tRNA|pre_miRNA|rRNA|miRNA/) && ($chromosome =~ m/^[[1-9]|10]$/)) { #specifies feature types to save, and assumes 10 chromosomes
	
		$id =~ m/gene:(.*?);/;
		my $symbol = $1;
		print $output "$chromosome\t$type\t$start\t$stop\t$symbol\n";
		
	} elsif (($type =~ m/mRNA/) && ($chromosome =~ m/^[[1-9]|10]$/)) {
	
		$id =~ m/transcript:(.*?);/;
		my $symbol = $1;
		$symbol =~ s/_T.*?$//;
		print $output "$chromosome\tgene\t$start\t$stop\t$symbol\n";
		
	}

}

close $input;
close $output;