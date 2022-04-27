#reanalysis and annotation scripts for GWAS results (will accept any properly formatted files, e.g. from GAPIT)
#useful for re-evaluating significance with new thresholds, or switching between user-defined, Bonferroni, and Benjamini-Hochberg versions and their combinations.

#configuration options start at line 34 for:
#	threshold options for selecting candidate significant markers. which of these options to use is defined at lines 202-204
#	method for assigning markers to genes-- specified LD window or n closest genes to marker

#STRUCTURE
#SD3/REANNOTATE/
#├── SD3-S1.pl #script to re-analyze GAPIT statistics files, generating new lists of significant markers
#├── input *_Statistics.txt file(s) containing GAPIT results 
#├── output *_Significance_Thresholds.txt file #significance thresholds 
#├── output *_SignificantHits.txt #significant markers
#├── output *_Arabidopsis.txt file #significant gene candidate matches to Arabidopsis
#├── output *_Rice.txt file #significant gene candidate matches to Rice
#└── COREFILES/
#	 ├── GOODMAN/1.GM.txt - 10.GM.txt #Goodman panel marker data in GAPIT numerical format (SD1-S2, SD1-S3)
#	 ├── NCRPIS/1.GM.txt - 10.GM.txt #NCRPIS panel marker data in GAPIT numerical format (SD1-S2, SD1-S3)
#	 ├── 3.2.1/1.GM.txt - 10.GM.txt #Goodman panel high-density marker data in GAPIT numerical format (SD1-S2, SD1-S3)
#	 ├── Rice.txt, Arabidopsis.txt #sequence similarity databases relating maize to rice and arabidopsis (SD1-S8, SD1-S9)
#	 ├── GeneNameList.gff3 #filtered gene info list (SD1-S5)
#	 ├── GenePositions.txt #list of gene positions (SD1-S6)
#	 ├── KnownGenes.txt #list of known genes in maize
#	 ├── Nearest.txt #atlas of nearest ten genes to each marker (SD1-S7)
#	 └── Manhattan.r #template for creating manhattan plot using R package manhattanly

use strict;
use warnings;
use File::Basename;
use List::Util qw(sum);
use POSIX qw(ceil);

my $version = 'GOODMAN'; #which genome to analyze - either GOODMAN, NCRPIS, or 3.2.1
my $bjh_alpha = 0.1; #Significance threshold for BJH FDR p-val
my $bonf_alpha = 0.05; #Significance threshold for Bonferroni alpha
my $p_annotate = 0.0001; #Significance threshold to annotate markers anyway

my $use_ld = 0; #whether to use number of close genes or LD window. $number_of_close_genes option will still function, e.g. if >1 gene within the LD window, so set appropriately.
my $number_of_close_genes = 10; #number of genes near each marker to annotate (max 10)
my $ld_window = 10000; #window in bp around each marker to find genes to annotate

my $type = fileparse_set_fstype("Unix"); #Allows spaces in file paths
my $dirname = dirname(__FILE__);
my $top = defined($dirname) ? $dirname : '.';

my @num_markers;
for my $i (1..10) { #Get numbers of markers on each chromosome

	open my $file, '<', "$top/COREFILES/$version/$i.GM.txt";
	while (<$file>) {}
	my $line = $.-1; #0-index takes care of -1 for empty line at end; another -1 for header
	push @num_markers, $line;
	close $file;

}

my @bonferroni;
my @pcalc;
for my $nm (@num_markers) {

	push @pcalc, $bonf_alpha/$nm; #Y-values for significance thresholds on plot
	my $x = -(log($bonf_alpha/$nm)/log(10));
	push @bonferroni, $x;
	
}

my $avg_p = sum(@pcalc)/@pcalc;
my $avg_Y = -(log($avg_p)/log(10)); #Average Bonferroni significance threshold across all 10 chromosomes


my %nonprot;
open my $gff, '<', "$top/COREFILES/GeneNameList.gff3";
while (<$gff>) {

	chomp $_;
	my @columns = split("\t", $_);
	if ($columns[1] ne 'gene') { $nonprot{$columns[4]} = $columns[1]; }

}

close $gff;

my %start;
my %stop;
if ($use_ld == 1) {

	open my $gff, '<', "$top/COREFILES/GeneNameList.gff3";
	while (<$gff>) {

		chomp $_;
		my @columns = split("\t", $_);
		
		if ($columns[1] eq 'gene') { 
		
			my $name = substr($columns[4], 0, -5);
			if ((!(exists($start{$name}))) || ($columns[2] < $start{$name})) { $start{$name} = $columns[2]; }
			if ((!(exists($stop{$name}))) || ($columns[3] > $stop{$name}))  { $stop{$name} = $columns[3]; }
		
		} else {
			
			if ((!(exists($start{$columns[4]}))) || ($columns[2] < $start{$columns[4]})) { $start{$columns[4]} = $columns[2]; }
			if ((!(exists($stop{$columns[4]}))) || ($columns[3] > $stop{$columns[4]}))  { $stop{$columns[4]} = $columns[3]; }
			
		}

	}

}

my %midpoints;
open my $mid, '<', "$top/COREFILES/GenePositions.txt";
while (<$mid>) {

	chomp $_;
	my @columns = split("\t", $_);
	$midpoints{$columns[2]} = ceil($columns[1]) ;

}

close $mid;


my %mappings;
open my $map, '<', "$top/COREFILES/Nearest.txt";
while (<$map>) {

	chomp $_;
	my @columns = split("\t", $_);
	if ($number_of_close_genes > 10) { $number_of_close_genes = 10; }
	$mappings{$columns[0]} = join("\t",@columns[1..$number_of_close_genes]);

}

close $map;

my %known;
open my $knownfile, '<', "$top/COREFILES/KnownGenes.txt";
my $header = <$knownfile>;
while (<$knownfile>) {

	chomp $_;
	my @columns = split("\t", $_);
	$known{$columns[0]} = join('/',@columns[1..2]);

}

close $knownfile;

my %rice;
open my $ricefile, '<', "$top/COREFILES/Rice.txt";
my $header1 = <$ricefile>;
while (<$ricefile>) {

	chomp $_;
	my @columns = split("\t", $_);
	$rice{$columns[0]} = join("\t",@columns[1..$#columns]);

}

close $ricefile;

my %arabidopsis;
open my $arabidopsisfile, '<', "$top/COREFILES/Arabidopsis.txt";
my $header2 = <$arabidopsisfile>;
while (<$arabidopsisfile>) {

	chomp $_;
	my @columns = split("\t", $_);
	$arabidopsis{$columns[0]} = join("\t",@columns[1..$#columns]);

}

close $arabidopsisfile;

my @Stats = glob("$top/*Statistics.txt");

for my $file (@Stats) {

$file =~ m/(.*?)_Statistics.txt/;
my $trait = $1;

open my $in, '<', "$top/$trait"."_Statistics.txt";
open my $out, '>', "$top/$trait"."_SignificantHits.txt";
my $head = (<$in>);

while (<$in>) { #Prints Significant Hits to output file

	chomp $_;
	my @names = split("\t", $_);
	my $SNP=$names[0];
	my $Chrom=$names[1];
	my $BP=$names[2];
	my $pval=$names[3];
	my $fdrpval=$names[8];
	
	#if (($fdrpval < $bjh_alpha) || ($pval < $pcalc[$Chrom-1]) || ($pval < $p_annotate)) { #must be below BJH, bonferroni, AND user-specified threshold
	if ($fdrpval < $bjh_alpha) { #must be below only BJH threshold
	#if ($pval < $pcalc[$Chrom-1]) { #must be below only user threshold
	
		my $tmpos = join(",",$SNP,$Chrom,$BP,$pval,$fdrpval);
		chomp $tmpos;
		print $out "$tmpos\n";
		
	}
	
}

close $out;
close $in;

#Find maize gene names from significant hit marker positions
open $in, '<', "$top/$trait"."_SignificantHits.txt";
open my $finalarabidopsis, '>', "$top/$trait"."_Arabidopsis.txt";
open my $finalrice, '>', "$top/$trait"."_Rice.txt";

print $finalarabidopsis "SNP\tChromosome\tLocus\tMLM p-value\tBJH p-value\tGene Midpoint\tGene Distance\tType\tMaize Gene\tArabidopsis Gene\tDescription\tBLAST Alignent Score\tE-value\n";
print $finalrice "SNP\tChromosome\tLocus\tMLM p-value\tBJH p-value\tGene Midpoint\tGene Distance\tType\tMaize Gene\tRice Gene\tDescription\tBLAST Alignent Score\tE-value\n";

while (<$in>) {

	chomp $_;
	my @gwaspos = split(',',$_);
	my $SNP=$gwaspos[0];
	my $chromosome=$gwaspos[1];
	my $posit=$gwaspos[2];
	my $pval=$gwaspos[3];
	my $fdrpval=$gwaspos[4];

	my @genelist = split("\t", $mappings{$SNP});
	for my $gene (@genelist) {
		
		if ($use_ld == 1) {

			my $startdistance = abs($start{$gene} - $posit);
			my $stopdistance = abs($stop{$gene} - $posit);
			unless (($startdistance <= $ld_window) || ($stopdistance <= $ld_window)) { next; }

		}
	
		my $midpoint = $midpoints{$gene};
		my $distance = $midpoint-$posit;
		my $type = $nonprot{$gene} // "gene";
		my $name = $known{$gene} // "none";
		my $rice = "none\tnone\tnone\tnone";
		my $arabidopsis = "none\tnone\tnone\tnone";
		
		if (exists($arabidopsis{$gene})) {
	
			$arabidopsis = $arabidopsis{$gene};

		}
		
		if (exists($rice{$gene})) {
		
			$rice = $rice{$gene};
		
		}
		
		print $finalarabidopsis "$SNP\t$chromosome\t$posit\t$pval\t$fdrpval\t$midpoint\t$distance\t$type\t$gene\t$name\t$arabidopsis\n";
		print $finalrice "$SNP\t$chromosome\t$posit\t$pval\t$fdrpval\t$midpoint\t$distance\t$type\t$gene\t$name\t$rice\n";
		
	}
	
	print $finalarabidopsis "\n";
	print $finalrice "\n";
	
}

close $in;
close $finalarabidopsis;
close $finalrice;

open $out, '>', "$top/$trait"."_Significance_Thresholds.txt";
print $out "Bonferroni Alpha: $bonf_alpha\nFDR-corrected p-value Threshold: $bjh_alpha\nAverage Bonferroni Threshold p-value: $avg_p\nAverage Bonferroni Threshold Y-value: $avg_Y\n";
close $out;

}