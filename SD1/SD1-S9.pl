#second of two scripts to generate sequence similarity databases between query and target proteomes, preceded by SD1-S9.pl
#requires the NCBI BLAST+ toolkit, BLAST databases formatted by said toolkit, and the gene names column from GeneNameList.gff3 (product of SD1-S5.pl) formatted as a standalone file

#target species specified line 31

#STRUCTURE
#SD1/
#├── SD1-S8.pl #first script to run to generate database
#├── SD1-S9.pl #second script to run to generate database
#├── GeneNameList.gff3 input file (product of SD1-S5.pl)
#├── SPECIES.txt #output file(s) for sequence similarity databases between query and target species
#└── SD1-Files/
#    └── BLAST/
#    	 ├── Various .exe and .dll files used in the NCBI BLAST+ toolkit
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

my @sp = ('Arabidopsis','Rice');
for my $species (@sp) {#Match Maize Query Genes to Rice/Arabidopsis BLAST Matches

	my %unique;
	open my $in, '<', "$top/Blasted$species.txt";
	open my $out, '>', "$top/$species.txt";
	print $out "Query\t$species Gene\tDescription\tBLAST Alignment Score\tE-Value\n";
	
	while (<$in>) {
		
			my @fields = split("\t",$_);
			my $gene = $fields[0];
		
			unless (exists($unique{$gene})) {
		
			$unique{$gene} = 1;
			$fields[1] =~ s/LOC_//;
			$fields[2] =~ s/%2C/,/; #Fix errors in BLAST descriptions from HTML-encoded characters
			$fields[2] =~ s/%26#64257%3B/fi/;
			$fields[2] =~ s/%26/&/;
			$fields[2] =~ s/%3B/;/;
			$fields[2] =~ s/\| (.*?)\|.*/$1/g;
			$fields[2] =~ s/protein\|//;
			my $description = $fields[2];
		
			my $outline = join("\t",$fields[0],$fields[1],$description,$fields[3],$fields[4]);
			chomp $outline;
			print $out "$outline\n"; 
		
		}
		
	}
		
	close $in;
	close $out;

}