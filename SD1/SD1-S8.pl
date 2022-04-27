#first of two scripts to generate sequence similarity databases between query and target proteomes, followed by SD1-S9.pl
#requires the NCBI BLAST+ toolkit, BLAST databases formatted by said toolkit, and the gene names column from GeneNameList.gff3 (product of SD1-S5.pl) formatted as a standalone file

#threshold for a valid match (i.e. value beyond which no match will be output to the final database) is an e-value of 1E-20 (specified in the '-evalue' parameters in lines 75/76 of SD1-S9.pl)
#if the best match value is above 1E-20, no match will be described. For more permissive matching, set the threshold higher in lines 75 and 76.

#query species specified line 51, and in FILES/BLAST/format_dbs.pl (when adding new databases)
#target species specified line 78, and in FILES/BLAST/format_dbs.pl (when adding new databases)
#name filtering options (specific to maize for this study) set on lines 58, 59

#STRUCTURE
#SD1/
#├── SD1-S8.pl #first script to run to generate database
#├── SD1-S9.pl #second script to run to generate database
#├── GeneNameList.gff3 input file (product of SD1-S5.pl)
#├── BlastedSPECIES.txt #output file(s) for sequence similarity databases between query and target species
#└── SD1-Files/
#    └── BLAST/
#    	 ├── .exe and .dll files used in the NCBI BLAST+ toolkit
#   	 ├── format_dbs.pl #run to use the NCBI BLAST+ toolkit to format new databases. Folder structure is DATABASES/SPECIES/SPECIES(.fasta)
#		 └── DATABASES/
#    		 └── folders for each species containing the FASTA-formatted proteome for that species, named identically to the folder and species within format_dbs.pl, this script, and SD1-S9.pl

use strict;
use warnings;
use Cwd qw(abs_path);
use File::Copy;
use File::Path;
use File::Basename;
use List::Util qw(sum);

my $type = fileparse_set_fstype("Unix"); #Allows spaces in file paths
my $dirname = dirname(__FILE__);
my $top = defined($dirname) ? $dirname : '.';

my @Hits;

open my $gff, '<', "$top/GeneNameList.gff3";
while (<$gff>) {
	
	chomp $_;
	my @GFFline = split("\t",$_);
	my $name=$GFFline[4];
	push @Hits, $name;

}

close $gff;

system("$top/SD1-Files/BLAST/blastdbcmd.exe -db $top/SD1-Files/BLAST/DATABASES/Maize/Formatted_Maize -dbtype prot -entry_batch $top/Names.txt -out $top/NamesRaw.txt");

open my $in, '<', "$top/NamesRaw.txt"; #Removes everything but the name from the ID line
open my $out, '>', "$top/NamesClean.txt";

while (<$in>) {

	my $linein = substr($_,0,1);
	my $newline = substr($_,1,14);

	if ($linein eq '>') {
	
		print $out '>'."$newline\n";
		
	} else {
	
		print $out $_;
		
	}
}
	
close $in;
close $out;

my @targetpecies = ("Maize","Rice");
for my $taxon (@targetspecies) {

	system("$top/SD1-Files/BLAST/blastp.exe -db $top/SD1-Files/BLAST/DATABASES/$taxon/Formatted_$taxon -query $top/NamesClean.txt -evalue 0.0000000000000000001 -out $top/Blasted$taxon.txt -outfmt \"6 qseqid sseqid stitle bitscore evalue\"");

}