#reanalysis and annotation scripts for GWAS results (will accept any properly formatted files, e.g. from GAPIT)
#useful for generating Manhattan plots with precisely placed vertical lines at marker or gene locations

#SD3/REANNOTATE/
#├── SD3-S2.pl #script to generate annotated Manhattan plots from GAPIT statistics files and list of genes
#├── input *_Statistics.txt file(s) #file containing GAPIT results 
#├── input *_SignificantHits.txt file(s) #file containing COMPILE-annotated significant results
#├── input *_Genes.txt file(s) #file containing list of genes to annotate
#├── output *_Manhattan.pdf file #annotated Manhattan plot
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
use Cwd qw(abs_path);
use File::Copy;
use File::Path;
use File::Basename;
use List::Util qw(sum);
use POSIX qw(ceil);

my $focus_mode = 0; #whether to analayze input files as regions of chromosomes, as SD2-S3 produces outputs
my $RPath = 'C:/"Program Files"/R/R-4.1.1/bin/x64'; #path to R installation, directories w/ spaces in double quotes
my $version = 'NCRPIS'; #which genome to analyze; either GOODMAN or NCRPIS
my $bjh_alpha = 0.1; #Significance threshold for BJH FDR p-val
my $bonf_alpha = 0.1; #Significance threshold for Bonferroni alpha

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

my %last_coords;
for my $i (1..10) {

	open my $fh, '<', "$top/COREFILES/$version/$i.GM.txt";
	my $lastline;
	$lastline = $_ while <$fh>;
	chomp $lastline;
	$last_coords{$i} = (split("\t",$lastline))[2];

}

my %chrom;
my %pos;

for my $i (1..10) {

	open my $fh, '<', "$top/COREFILES/$version/$i.GM.txt";
	my $header = <$fh>;
	
	while (<$fh>) {
	
		chomp $_;
		my @columns = split("\t",$_);
		$chrom{$columns[0]} = $columns[1];
		$pos{$columns[0]} = $columns[2];
	
	}

}

my %firsts;
my $counter = 1;

for my $i (1..10) {

$firsts{$i} = $counter;
my $foo = $last_coords{$i};
$counter += $last_coords{$i};

}

my @bonferroni;
my @pcalc;
for my $nm (@num_markers) {

	push @pcalc, $bonf_alpha/$nm; #Y-values for significance thresholds on plot
	my $x = -(log($bonf_alpha/$nm)/log(10));
	push @bonferroni, $x;
	
}

my $b1 = $bonferroni[0]; #Because R doesn't accept the array value; I know it looks bad
my $b2 = $bonferroni[1];
my $b3 = $bonferroni[2];
my $b4 = $bonferroni[3];
my $b5 = $bonferroni[4];
my $b6 = $bonferroni[5];
my $b7 = $bonferroni[6];
my $b8 = $bonferroni[7];
my $b9 = $bonferroni[8];
my $b10 = $bonferroni[9];

my $avg_p = sum(@pcalc)/@pcalc;
my $avg_Y = -(log($avg_p)/log(10)); #Average Bonferroni significance threshold across all 10 chromosomes

my %midpoints;
my %chromosomes;
open my $mid, '<', "$top/COREFILES/GenePositions.txt";
while (<$mid>) {

	chomp $_;
	my @columns = split("\t", $_);
	$midpoints{$columns[2]} = ceil($columns[1]) ;
	$chromosomes{$columns[2]} = $columns[0];

}

close $mid;

my @Stats = glob("$top/*Statistics.txt");

for my $file (@Stats) {

$file =~ m/\.\/(.*?)_Statistics.txt/;
my $trait = $1;

open my $in,  '<', "$top/$trait"."_Statistics.txt";
open my $out, '>', "$top/$trait"."_StatisticsTMP.txt";
print $out "SNP\tChromosome\tPosition\tP.value\tx\ttaxa\tx\ty\tFDR.p.value\n";

my @coords;
my $current_chrom;
if ($focus_mode = 0) {

	while (<$in>) {

		print $out $_;
	
	}
	
} else {
	
	my $i = 0;
	while (<$in>) {

		push @coords, (split("\t", $_))[2];
		print $out $_;
	
		if ($i == 0) { 
	
			$i++; 
			$current_chrom = (split("\t", $_))[1]
		
		}
	
	}

}

close $in;
close $out;

my @sorted_coords;
my $pos1 $pos2;
if ($focus_mode = 1) { 

	my @sorted_coords = sort { $a <=> $b } @coords; }
	my $pos1 = shift(@sorted_coords)/1000000;
	my $pos2 = pop(@sorted_coords)/1000000;
	
	open my $manhattan, '<', "$top/COREFILES/focus_manhattan.r";
	my $file_content = do { local $/; <$manhattan> };
	open my $local_manhattan, '>', "$top/COREFILES/TMP_manhattan.r";
	$file_content =~ s/insertxmin/xmin = $pos1/;
	$file_content =~ s/insertxmax/xmax = $pos2/;
	print $local_manhattan $file_content;
	close $manhattan;
	close $local_manhattan;
	
}
	
my @hits;
open my $tmp,  '<', "$top/$trait"."_SignificantHits.txt";

while (<$tmp>) {

	push @hits, (split(',', $_))[0]; 
	
}

close $tmp;

my @genes;
if (-e "$top/$trait"."_genes.txt") {
	open my $genelist, '<', "$top/$trait"."_genes.txt";

	while (<$genelist>) {

		chomp $_;
		push @genes, $_; 
	
	}

	close $genelist;

}

my @markers;
if (-e "$top/$trait"."_markers.txt") {
	open my $markerlist, '<', "$top/$trait"."_markers.txt";

	while (<$markerlist>) {

		chomp $_;
		push @markers, $_; 
	
	}

	close $markerlist;

}

my %strings;
my $total_string;

my $color = qq("red");

if ($focus_mode = 0) {

foreach my $gene (@genes) {

	my $rawposition = $midpoints{$gene};
	my $position = $rawposition/1000000;
	my $chromosome = $chromosomes{$gene};
	my $total_position = $rawposition + $firsts{$chromosome};
	
	if (exists($strings{$chromosome})) { $strings{$chromosome} = qq($strings{$chromosome},"$position"); }
	else { $strings{$chromosome} = qq("$position"); }
	
	if (defined($total_string)) { $total_string = qq($total_string,"$total_position"); }
	else { $total_string = qq("$total_position"); }
	
}

foreach my $marker (@markers) {

	my $rawposition = $pos{$marker};
	my $position = $rawposition/1000000;
	my $chromosome = $chrom{$marker};
	my $total_position = $rawposition + $firsts{$chromosome};
	
	if (exists($strings{$chromosome})) { $strings{$chromosome} = qq($strings{$chromosome},"$position"); }
	else { $strings{$chromosome} = qq("$position"); }
	
	if (defined($total_string)) { $total_string = qq($total_string,"$total_position"); }
	else { $total_string = qq("$total_position"); }
	
}

} else {
	
foreach my $gene (@genes) {

	my $rawposition = $midpoints{$gene};
	my $position = $rawposition/1000000;
	my $chromosome = $chromosomes{$gene};
	
	if (defined($total_string)) { $total_string = qq($total_string,"$position"); }
	else { $total_string = qq("$position"); }
	
}
	
foreach my $marker (@markers) {

	my $rawposition = $midpoints{$marker};
	my $position = $rawposition/1000000;
	my $chromosome = $chromosomes{$marker};
	
	if (defined($total_string)) { $total_string = qq($total_string,"$position"); }
	else { $total_string = qq("$position"); }
	
}

}

my $string1 = '';
my $string2 = '';
my $string3 = '';
my $string4 = '';
my $string5 = '';
my $string6 = '';
my $string7 = '';
my $string8 = '';
my $string9 = '';
my $string10 = '';

if (defined($total_string)) { $total_string = "$total_string"; }
else { $total_string = ''; }

if (exists($strings{1})) { $string1 = "$strings{1}"; }
else { $string1 = ''; }

if (exists($strings{2})) { $string2 = "$strings{2}"; }
else { $string2 = ''; }

if (exists($strings{3})) { $string3 = "$strings{3}"; }
else { $string3 = ''; }

if (exists($strings{4})) { $string4 = "$strings{4}"; }
else { $string4 = ''; }

if (exists($strings{5})) { $string5 = "$strings{5}"; }
else { $string5 = ''; }

if (exists($strings{6})) { $string6 = "$strings{6}"; }
else { $string6 = ''; }

if (exists($strings{7})) { $string7 = "$strings{7}"; }
else { $string7 = ''; }

if (exists($strings{8})) { $string8 = "$strings{8}"; }
else { $string8 = ''; }

if (exists($strings{9})) { $string9 = "$strings{9}"; }
else { $string9 = ''; }

if (exists($strings{10})) { $string10 = "$strings{10}"; }
else { $string10 = ''; }

my $R_HITS = join('","', @hits);

my $core = dirname(abs_path($0));
open my $rfile, '>', "$top/R.txt";

if ($focus_mode = 0) { 

print $rfile qq`capture.output()
setwd("$core")
name <- "$trait"
library(compiler)
library(gplots)
library(genetics)
library(EMMREML)
library("scatterplot3d")
library(MASS)
source("$core/COREFILES/manhattan.r")
library(manhattanly)
library(plotly)
GWAS.Results <- read.table(paste("$core/$trait", "_StatisticsTMP.txt", sep = ""), head=TRUE)
SNPsOfInterest <- c("$R_HITS")
manhattan(GWAS.Results, genomewideline=$avg_Y, chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black", "gold3"), main = "Manhattan Plot for $trait", highlight = SNPsOfInterest, ablines=c($total_string))
manhattan(subset(GWAS.Results, Chromosome==1), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b1, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest, ablines=c($string1))
manhattan(subset(GWAS.Results, Chromosome==2), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b2, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest, ablines=c($string2))
manhattan(subset(GWAS.Results, Chromosome==3), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b3, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest, ablines=c($string3))
manhattan(subset(GWAS.Results, Chromosome==4), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b4, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest, ablines=c($string4))
manhattan(subset(GWAS.Results, Chromosome==5), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b5, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest, ablines=c($string5))
manhattan(subset(GWAS.Results, Chromosome==6), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b6, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest, ablines=c($string6))
manhattan(subset(GWAS.Results, Chromosome==7), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b7, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest, ablines=c($string7))
manhattan(subset(GWAS.Results, Chromosome==8), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b8, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest, ablines=c($string8))
manhattan(subset(GWAS.Results, Chromosome==9), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b9, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest, ablines=c($string9))
manhattan(subset(GWAS.Results, Chromosome==10), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b10, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest, ablines=c($string10))
q()
`;

} else {
	
print $rfile qq`capture.output()
setwd("$core")
name <- "$trait"
library(compiler)
library(gplots)
library(genetics)
library(EMMREML)
library("scatterplot3d")
library(MASS)
source("$top/COREFILES/TMP_manhattan.r")
library(manhattanly)
library(plotly)
GWAS.Results <- read.table(paste("$core/$trait", "_StatisticsTMP.txt", sep = ""), head=TRUE)
SNPsOfInterest <- c("$R_HITS")
manhattan(subset(GWAS.Results, Chromosome==$current_chrom), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$bon, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest, ablines=c($total_string))
q()
`;
	
}

close $rfile;
my $command = "$RPath/Rscript $top/R.txt";
system("$command");

unlink("$top/$trait"."_StatisticsTMP.txt");
unlink("$top/R.txt");
unlink("$top/COREFILES/TMP_manhattan.r");
move("$top/Rplots.pdf", "$top/$trait"."_Manhattan.pdf");

}