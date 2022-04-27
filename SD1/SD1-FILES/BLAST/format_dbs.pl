use strict;
use warnings;

my @species = ("X","Y","Z");

for my $taxon (@species) {
	
	system("makeblastdb.exe -in DATABASES/$taxon/$taxon.fasta -dbtype prot -parse_seqids -title Formatted_$taxon -out DATABASES/$taxon/Formatted_$taxon");
	
}