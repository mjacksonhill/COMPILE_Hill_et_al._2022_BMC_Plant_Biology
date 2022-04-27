#configuration options start at line 38 for:
#	threshold options for selecting candidate significant markers. which of these options to use is defined at lines 484-485 and 571-573
#	method for assigning markers to genes-- specified LD window or n closest genes to marker

#more in-depth customization:
#	GAPIT model defined at lines 343, 381, together with surrounding script. For example, enable compression by setting the group.from parameter to 0, and/or by adding a group.to parameter to adjust the groupings.

#STRUCTURE
#SD2/GOODMAN_2.7/
#├── SD2-S1.pl #script to execute COMPILE
#├── input *.txt file(s) containing phenotype data
#├── RUNS/
#│	 ├── Example.txt #contains example input data in correct format
#│	 └── RESULTS_FOLDERS #one per input data file 
#│		 ├── INTERMEDIATES/ #raw GAPIT output files, including some not used by COMPILE (e.g. reports of marker effect size)
#│		 ├── output *_DATA.txt #copy of original input phenotype data
#│		 ├── output *_Manhattan.pdf #Manhattan plots of results
#│		 ├── output *_QQ.pdf #quantile-quantile plot of results
#│		 ├── output *_Arabidopsis.txt #candidate gene report matching maize to arabidopsis genes
#│		 └── output *_Rice.txt #candidate gene report matching maize to rice genes
#└── COREFILES/
#	 ├── GAPIT/ #contains files used by GAPIT to execute GWAS
#	 └── GENOME/
#		 ├── 1.GD.txt - 10.GD.txt, 1.GM.txt - 10.GM.txt #marker data in GAPIT numerical format (SD1-S2, SD1-S3)
#		 ├── 1k.txt - 10k.txt #kinship matrix data (SD1-S4)
#		 ├── GeneNameList.gff3 #filtered gene info list (SD1-S5)
#		 ├── GenePositions.txt #list of gene positions (SD1-S6)
#		 ├── KnownGenes.txt #list of known genes in maize
#		 ├── taxa.txt #list of taxa
#		 ├── Nearest.txt #atlas of nearest ten genes to each marker (SD1-S7)
#		 └── Rice.txt, Arabidopsis.txt #sequence similarity databases relating maize to rice and arabidopsis (SD1-S8, SD1-S9)

use strict;
use warnings;
use Cwd qw(abs_path);
use File::Copy;
use File::Path;
use File::Basename;
use List::Util qw(sum);
use POSIX qw(ceil);

my $debug = 0; #If 1, doesn't delete temporary files or directories
my $do_serial = 0; #If 1, runs chromosomes in order 1-10 instead of in parallel
my $bjh_alpha = 0.1; #Significance threshold for BJH FDR p-val
my $bonf_alpha = 0.1; #Significance threshold for Bonferroni alpha
my $p_annotate = 0.0001; #Significance threshold to annotate markers anyway

my $use_ld = 1; #whether to use number of close genes or LD window. $number_of_close_genes option will still function, e.g. if >1 gene within the LD window, so set appropriately.
my $number_of_close_genes = 10; #number of genes near each marker to annotate (max 10)
my $ld_window = 10000; #window in bp around each marker to find genes to annotate

my $type = fileparse_set_fstype("Unix"); #Allows spaces in file paths
my $dirname = dirname(__FILE__);
my $top = defined($dirname) ? $dirname : '.';

my @num_markers;
for my $i (1..10) { #Get numbers of markers on each chromosome

	open my $file, '<', "$top/COREFILES/GENOME/$i.GM.txt";
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

print "Initializing, please wait a moment...\n";

my %nonprot;
open my $gff, '<', "$top/COREFILES/GENOME/GeneNameList.gff3";
while (<$gff>) {

	chomp $_;
	my @columns = split("\t", $_);
	if ($columns[1] ne 'gene') { $nonprot{$columns[4]} = $columns[1]; }

}

close $gff;

my %start;
my %stop;
if ($use_ld == 1) {

	open my $gff, '<', "$top/COREFILES/GENOME/GeneNameList.gff3";
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
open my $mid, '<', "$top/COREFILES/GENOME/GenePositions.txt";
while (<$mid>) {

	chomp $_;
	my @columns = split("\t", $_);
	$midpoints{$columns[2]} = ceil($columns[1]) ;

}

close $mid;

my %mappings;
open my $map, '<', "$top/COREFILES/GENOME/Nearest.txt";
while (<$map>) {

	chomp $_;
	my @columns = split("\t", $_);
	if ($number_of_close_genes > 10) { $number_of_close_genes = 10; }
	$mappings{$columns[0]} = join("\t",@columns[1..$number_of_close_genes]);

}

close $map;

my %known;
open my $knownfile, '<', "$top/COREFILES/Genome/KnownGenes.txt";
my $header1 = <$knownfile>;
while (<$knownfile>) {

	chomp $_;
	my @columns = split("\t", $_);
	$known{$columns[0]} = join('/',@columns[1..2]);

}

close $knownfile;

my %rice;
open my $ricefile, '<', "$top/COREFILES/Genome/Rice.txt";
my $header = <$ricefile>;
while (<$ricefile>) {

	chomp $_;
	my @columns = split("\t", $_);
	$rice{$columns[0]} = join("\t",@columns[1..$#columns]);

}

close $ricefile;

my %arabidopsis;
open my $arabidopsisfile, '<', "$top/COREFILES/Genome/Arabidopsis.txt";
my $header2 = <$arabidopsisfile>;
while (<$arabidopsisfile>) {

	chomp $_;
	my @columns = split("\t", $_);
	$arabidopsis{$columns[0]} = join("\t",@columns[1..$#columns]);

}

close $arabidopsisfile;

open my $RPath, '<', "$top/COREFILES/R_Installation.txt";
my $R = <$RPath>;
close $RPath;

print "Please place copies of your phenotype data in .txt format in the base Pipeline folder.\nThe Pipeline will run them sequentially.\nData files will be renamed and moved to the run folder in /RUNS/.\n";
print "Please ensure the phenotype data is in the format in /RUNS/Template.txt.\nPress any key to continue.\n>";
my $GO = <>;

my @datafiles = glob("$top/*.txt");
for my $file (@datafiles) { #Executes the Pipeline for every set of trait data in the folder
	
opendir my $dh, "$top/RUNS"; #Gets the name of every existing directory in RUNS
my @existing_runs = grep {-d "$top/RUNS/$_" && ! /^\.{1,2}$/} readdir($dh);
closedir $dh;

open my $traitfile, '<', "$file"; #Gets the first line of the input file
my $firstline = <$traitfile>;
close $traitfile;

$firstline =~ /\t(.*)$/; #Gets the trait name from the first line
my $trait = $1;
chomp $trait;
my $traitund = "$trait"."__";

my $numfix = 1; #Finds how many runs for this trait already exist in RUNS and add 1 to the suffix
for my $run (@existing_runs) {

	my $suffx = (split('__',$run))[1];
	if ($run =~ m/$traitund/) { 
	
		$numfix = $suffx+1; 
	
	}
	
}

my $RunName = join('__', $trait,$numfix);

mkdir "$top/RUNS/$RunName", 0755; #Creates and populates the run directory
mkdir "$top/RUNS/$RunName/TMP", 0755;
my $RunDataName = join('_', $trait,'DATA');

open my $translate, '<', "$file";
open my $trans_out, '>', "$top/RUNS/$RunName/$RunDataName.txt";
while (<$translate>) {

	#$_ =~ s///; #In case some global translation of input taxa names is needed in future
	print $trans_out $_;
	
}

close $trans_out;

open my $names, '<', "$top/RUNS/$RunName/$RunDataName.txt"; #Slurp in lines from data file
my @names = <$names>;
close $names;

my @exists; #Create an array of taxa with data
for (@names) {

	if ($_ !~ m/NaN/) { push @exists, (split("\t", $_))[0]; }
	
}

chomp @exists;

open my $reference, '<', "$top/COREFILES/GENOME/taxa.txt"; #Nab the taxa list
my @columns = split("\t", <$reference>); #Get an array of taxa names
close $reference;
chomp @columns;

my @keep; #Array of taxa names to keep

for my $j (0..$#columns) { #Pushes only taxa for which there is data to @keep by looping through @exists for each value in @columns

	my $flg = 0;
	for my $k (0..$#exists) {
	
		if ($exists[$k] eq $columns[$j]) { $flg = 1; } #If the taxon is present, set the flag to 1 so it gets pushed to @keep
	
	}
	if ($flg == 1) { push @keep, $columns[$j]; }

}

my $taxa = scalar(@keep); #Number of taxa is equal to the number of kept columns
my %to_keep = map { $_ => 1 } @keep;

print "Scaling Marker Sets ($taxa Individuals)...\n";
for my $i (1..10) {

	die "could not fork" unless defined(my $prepid = fork); #Parallel processing of genotype files
	unless ($prepid) {
	
		copy("$top/COREFILES/GENOME/$i.GM.txt", "$top/RUNS/$RunName/TMP/$i.GM.txt");
		open my $in, '<', "$top/COREFILES/GENOME/$i.GD.txt";
		open my $out, '>', "$top/RUNS/$RunName/TMP/$i.GD.txt";
		my $header = <$in>;
		print $out $header;
	
		while (<$in>) {
		
			my $taxon = (split("\t", $_))[0];
			if(exists($to_keep{$taxon})) { print $out $_; }	
		
		}
	
		close $in;
		close $out;
		exit;
		
	}

}

while (1) {

  my $child = waitpid(-1, 0);
  last if $child == -1;

}

if ($do_serial == 1) {

print "Executing Chromosomal MLM via GAPIT for $trait serially...\n";
for my $i (1..10) { #Simultaneously executes GAPIT for all 10 chromosomes
	
	my $k = "$i"."k";
	my $I = $i;
	mkdir "$top/RUNS/$RunName/TMP/$i", 0755;
	open my $rfile, '>', "$top/RUNS/$RunName/TMP/$i/R.txt";
	my $core = dirname(abs_path($0));
	
	print $rfile qq`library(compiler)
source("$core/COREFILES/GAPIT/gapit_functions.txt")
source("$core/COREFILES/GAPIT/emma.txt")
source("$core/COREFILES/GAPIT/mt.txt")
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(scatterplot3d)
myGD <- read.table("$core/RUNS/$RunName/TMP/$i.GD.txt", head=TRUE)
myGM <- read.table("$core/RUNS/$RunName/TMP/$i.GM.txt", head=TRUE)
myKI <- read.table("$core/COREFILES/GENOME/$k.txt", head=FALSE)
myY <- read.table("$core/RUNS/$RunName/$RunDataName.txt", head=TRUE)
setwd("$top/RUNS/$RunName/TMP/$i")
myGAPIT <- GAPIT(Y=myY, GD=myGD,GM=myGM, KI=myKI, group.from=282, group.to=282, Geno.View.output=FALSE)
q()
`;	
	close $rfile;
	my $command = "$R/Rscript $top/RUNS/$RunName/TMP/$i/R.txt > nul 2>\&1";
	system("$command");
	print "Chromosome $i complete!\n";

}

} else {

print "Executing Chromosomal MLM via GAPIT for $trait in parallel...\n";
for my $i (1..10) { #Simultaneously executes GAPIT for all 10 chromosomes

	die "could not fork" unless defined(my $pid = fork);
	unless ($pid) {
	
		my $k = "$i"."k";
		my $I = $i;
		mkdir "$top/RUNS/$RunName/TMP/$i", 0755;
		open my $rfile, '>', "$top/RUNS/$RunName/TMP/$i/R.txt";
		my $core = dirname(abs_path($0));
	
		print $rfile qq`library(compiler)
source("$core/COREFILES/GAPIT/gapit_functions.txt")
source("$core/COREFILES/GAPIT/emma.txt")
source("$core/COREFILES/GAPIT/mt.txt")
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(scatterplot3d)
myGD <- read.table("$core/RUNS/$RunName/TMP/$i.GD.txt", head=TRUE)
myGM <- read.table("$core/RUNS/$RunName/TMP/$i.GM.txt", head=TRUE)
myKI <- read.table("$core/COREFILES/GENOME/$k.txt", head=FALSE)
myY <- read.table("$core/RUNS/$RunName/$RunDataName.txt", head=TRUE)
setwd("$core/RUNS/$RunName/TMP/$i")
myGAPIT <- GAPIT(Y=myY, GD=myGD,GM=myGM, KI=myKI, group.from=282, group.to=282, Geno.View.output=FALSE)
q()
`;	
		close $rfile;
		my $command = "$R/Rscript $top/RUNS/$RunName/TMP/$i/R.txt > nul 2>\&1";
		system("$command");
		print "Chromosome $i complete!\n";
		exit;
		
	}

}

while (1) {

  my $child = waitpid(-1, 0);
  last if $child == -1;

}

}

print "MLM Complete! Compiling results...\n"; #Compiles GWAS Results

open my $rfile, '>', "$top/RUNS/$RunName/TMP/R1.txt";
my $core = dirname(abs_path($0));

print $rfile qq`library(compiler)
source("$core/COREFILES/GAPIT/gapit_functions.txt")
source("$core/COREFILES/GAPIT/mt.txt")
library(gplots)
library(genetics)
library(EMMREML)
library("scatterplot3d")
library(MASS)
source("$core/COREFILES/GAPIT/manhattan.txt")
library(manhattanly)
library(plotly)
mydataPath.Results.1="$core/RUNS/$RunName/TMP/1/"
mydataPath.Results.2="$core/RUNS/$RunName/TMP/2/"
mydataPath.Results.3="$core/RUNS/$RunName/TMP/3/"
mydataPath.Results.4="$core/RUNS/$RunName/TMP/4/"
mydataPath.Results.5="$core/RUNS/$RunName/TMP/5/"
mydataPath.Results.6="$core/RUNS/$RunName/TMP/6/"
mydataPath.Results.7="$core/RUNS/$RunName/TMP/7/"
mydataPath.Results.8="$core/RUNS/$RunName/TMP/8/"
mydataPath.Results.9="$core/RUNS/$RunName/TMP/9/"
mydataPath.Results.10="$core/RUNS/$RunName/TMP/10/"
name <- "$trait"
GWAS.Results.1 <- read.csv(paste(mydataPath.Results.1,"GAPIT..",name,".GWAS.Results.csv",sep=""), head=TRUE)
GWAS.Results.2 <- read.csv(paste(mydataPath.Results.2,"GAPIT..",name,".GWAS.Results.csv",sep=""), head=TRUE)
GWAS.Results.3 <- read.csv(paste(mydataPath.Results.3,"GAPIT..",name,".GWAS.Results.csv",sep=""), head=TRUE)
GWAS.Results.4 <- read.csv(paste(mydataPath.Results.4,"GAPIT..",name,".GWAS.Results.csv",sep=""), head=TRUE)
GWAS.Results.5 <- read.csv(paste(mydataPath.Results.5,"GAPIT..",name,".GWAS.Results.csv",sep=""), head=TRUE)
GWAS.Results.6 <- read.csv(paste(mydataPath.Results.6,"GAPIT..",name,".GWAS.Results.csv",sep=""), head=TRUE)
GWAS.Results.7 <- read.csv(paste(mydataPath.Results.7,"GAPIT..",name,".GWAS.Results.csv",sep=""), head=TRUE)
GWAS.Results.8 <- read.csv(paste(mydataPath.Results.8,"GAPIT..",name,".GWAS.Results.csv",sep=""), head=TRUE)
GWAS.Results.9 <- read.csv(paste(mydataPath.Results.9,"GAPIT..",name,".GWAS.Results.csv",sep=""), head=TRUE)
GWAS.Results.10 <- read.csv(paste(mydataPath.Results.10,"GAPIT..",name,".GWAS.Results.csv",sep=""), head=TRUE)
Effect.Estimates.1 <- read.csv(paste(mydataPath.Results.1,"GAPIT..",name,".Allelic_Effect_Estimates.csv",sep=""), head=TRUE)
Effect.Estimates.2 <- read.csv(paste(mydataPath.Results.2,"GAPIT..",name,".Allelic_Effect_Estimates.csv",sep=""), head=TRUE)
Effect.Estimates.3 <- read.csv(paste(mydataPath.Results.3,"GAPIT..",name,".Allelic_Effect_Estimates.csv",sep=""), head=TRUE)
Effect.Estimates.4 <- read.csv(paste(mydataPath.Results.4,"GAPIT..",name,".Allelic_Effect_Estimates.csv",sep=""), head=TRUE)
Effect.Estimates.5 <- read.csv(paste(mydataPath.Results.5,"GAPIT..",name,".Allelic_Effect_Estimates.csv",sep=""), head=TRUE)
Effect.Estimates.6 <- read.csv(paste(mydataPath.Results.6,"GAPIT..",name,".Allelic_Effect_Estimates.csv",sep=""), head=TRUE)
Effect.Estimates.7 <- read.csv(paste(mydataPath.Results.7,"GAPIT..",name,".Allelic_Effect_Estimates.csv",sep=""), head=TRUE)
Effect.Estimates.8 <- read.csv(paste(mydataPath.Results.8,"GAPIT..",name,".Allelic_Effect_Estimates.csv",sep=""), head=TRUE)
Effect.Estimates.9 <- read.csv(paste(mydataPath.Results.9,"GAPIT..",name,".Allelic_Effect_Estimates.csv",sep=""), head=TRUE)
Effect.Estimates.10 <- read.csv(paste(mydataPath.Results.10,"GAPIT..",name,".Allelic_Effect_Estimates.csv",sep=""), head=TRUE)
GWAS.Results <- rbind(GWAS.Results.1, GWAS.Results.2, GWAS.Results.3, GWAS.Results.4, GWAS.Results.5, GWAS.Results.6, GWAS.Results.7, GWAS.Results.8, GWAS.Results.9, GWAS.Results.10)
Effect.Estimates <- rbind(Effect.Estimates.1, Effect.Estimates.2, Effect.Estimates.3, Effect.Estimates.4, Effect.Estimates.5, Effect.Estimates.6, Effect.Estimates.7, Effect.Estimates.8, Effect.Estimates.9, Effect.Estimates.10)
GWAS.Results <- GWAS.Results[,-ncol(GWAS.Results)]
setwd("$core/RUNS/$RunName")
Conduct.FDR <- GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(PWI=GWAS.Results, FDR.Rate=$bjh_alpha, FDR.Procedure="BH")
GWAS.Results.FDR <- Conduct.FDR\$PWIP
write.table(GWAS.Results.FDR, paste("GAPIT.", name, ".GWAS.Results.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(Effect.Estimates, paste("GAPIT.", name, ".Allelic_Effect_Estimates.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
`;

close $rfile;
my $command = "$R/Rscript $top/RUNS/$RunName/TMP/R1.txt > nul 2>\&1";
system("$command");

move("$top/RUNS/$RunName/GAPIT.$trait.Allelic_Effect_Estimates.txt", "$top/RUNS/$RunName/$trait"."_Effects.txt");
move("$top/RUNS/$RunName/GAPIT.$trait.GWAS.Results.txt", "$top/RUNS/$RunName/$trait"."_StatisticsTMP.txt");

open my $in, '<', "$top/RUNS/$RunName/$trait"."_StatisticsTMP.txt"; #Creates a new statistics file sorted by corrected FDR P-value
my $discard = <$in>;
my @lines = <$in>;
chomp @lines;
close $in;

my @ordered_lines = sort { (split("\t", $a))[1] <=> (split("\t", $b))[1] || (split("\t", $a))[2] <=> (split("\t", $b))[2] } @lines;

open my $out, '>', "$top/RUNS/$RunName/$trait"."_Statistics.txt";
for (@ordered_lines) { print $out "$_\n"; }
close $out;

my @BJH; #Pushes to @BJH FDR P-values below $bjh_alpha, these hits are colored in the final Manhattan plot
open $in,  '<',  "$top/RUNS/$RunName/$trait"."_Statistics.txt";

while (<$in>) {

	#if ((split("\t", $_))[8] < $bjh_alpha) { #bjh-corrected p-value below BJH threshold
	if ((split("\t", $_))[3] < $pcalc[(split("\t", $_))[1]-1]) { #uncorrected p-value below Bonferroni threshold
	
		push @BJH, (split("\t", $_))[0]; 
		
	}
	
}

close $in;
my $R_HITS = join('","', @BJH);

#Assimilate the final interactive Manhattan plots
open $rfile, '>', "$top/RUNS/$RunName/TMP/R2.txt";
print $rfile qq`capture.output()
setwd("$core/RUNS/$RunName")
name <- "$trait"
library(compiler)
source("$core/COREFILES/GAPIT/gapit_functions.txt")
source("$core/COREFILES/GAPIT/mt.txt")
library(gplots)
library(genetics)
library(EMMREML)
library("scatterplot3d")
library(MASS)
source("$core/COREFILES/GAPIT/manhattan.txt")
library(manhattanly)
library(plotly)
GWAS.Results <- read.table(paste("$core/RUNS/$RunName/$trait", "_StatisticsTMP.txt", sep = ""), head=TRUE)
SNPsOfInterest <- c("$R_HITS")
GAPIT.QQ(P.values = GWAS.Results[,4], name.of.trait = name,DPP=50000)
manhattan(GWAS.Results, genomewideline=$avg_Y, chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black", "gold3"), main = "Manhattan Plot for $trait", highlight = SNPsOfInterest)
manhattan(subset(GWAS.Results, Chromosome==1), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b1, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest)
manhattan(subset(GWAS.Results, Chromosome==2), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b2, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest)
manhattan(subset(GWAS.Results, Chromosome==3), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b3, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest)
manhattan(subset(GWAS.Results, Chromosome==4), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b4, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest)
manhattan(subset(GWAS.Results, Chromosome==5), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b5, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest)
manhattan(subset(GWAS.Results, Chromosome==6), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b6, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest)
manhattan(subset(GWAS.Results, Chromosome==7), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b7, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest)
manhattan(subset(GWAS.Results, Chromosome==8), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b8, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest)
manhattan(subset(GWAS.Results, Chromosome==9), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b9, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest)
manhattan(subset(GWAS.Results, Chromosome==10), chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black"), genomewideline=$b10, main = "Manhattan Plot for $trait", highlight = SNPsOfInterest)
GR <- manhattanr(GWAS.Results, chr = "Chromosome", bp = "Position", p = "P.value")
htmlwidgets::saveWidget(manhattanly(subset(GR[["data"]], CHR %in% 1), genomewideline=$b1, title = "Manhattan Plot for $trait", col=c("black")), "Chr_1.html")
htmlwidgets::saveWidget(manhattanly(subset(GR[["data"]], CHR %in% 2), genomewideline=$b2, title = "Manhattan Plot for $trait", col=c("black")), "Chr_2.html")
htmlwidgets::saveWidget(manhattanly(subset(GR[["data"]], CHR %in% 3), genomewideline=$b3, title = "Manhattan Plot for $trait", col=c("black")), "Chr_3.html")
htmlwidgets::saveWidget(manhattanly(subset(GR[["data"]], CHR %in% 4), genomewideline=$b4, title = "Manhattan Plot for $trait", col=c("black")), "Chr_4.html")
htmlwidgets::saveWidget(manhattanly(subset(GR[["data"]], CHR %in% 5), genomewideline=$b5, title = "Manhattan Plot for $trait", col=c("black")), "Chr_5.html")
htmlwidgets::saveWidget(manhattanly(subset(GR[["data"]], CHR %in% 6), genomewideline=$b6, title = "Manhattan Plot for $trait", col=c("black")), "Chr_6.html")
htmlwidgets::saveWidget(manhattanly(subset(GR[["data"]], CHR %in% 7), genomewideline=$b7, title = "Manhattan Plot for $trait", col=c("black")), "Chr_7.html")
htmlwidgets::saveWidget(manhattanly(subset(GR[["data"]], CHR %in% 8), genomewideline=$b8, title = "Manhattan Plot for $trait", col=c("black")), "Chr_8.html")
htmlwidgets::saveWidget(manhattanly(subset(GR[["data"]], CHR %in% 9), genomewideline=$b9, title = "Manhattan Plot for $trait", col=c("black")), "Chr_9.html")
htmlwidgets::saveWidget(manhattanly(subset(GR[["data"]], CHR %in% 10), genomewideline=$b10, title = "Manhattan Plot for $trait", col=c("black")), "Chr_10.html")
rm(GWAS.Results)
rm(Effect.Estimates)
`;

close $rfile;
$command = "$R/Rscript $top/RUNS/$RunName/TMP/R2.txt > nul 2>\&1";
system("$command");

print "Results Compiled! Cleaning...\n"; #Finish sorting out GAPIT outputs
move("$top/RUNS/$RunName/Rplots.pdf", "$top/RUNS/$RunName/$trait"."_Manhattan.pdf");
move("$top/RUNS/$RunName/GAPIT.$trait.QQ-Plot.pdf", "$top/RUNS/$RunName/$trait"."_QQ.pdf");

unless ($debug == 1) { rmtree("$top/RUNS/$RunName/TMP"); }
print "Results Cleaned! Beginning Interpretation...\n";

my @Hi_Mom;

open $out, '>', "$top/RUNS/$RunName/$trait"."_SignificantHits";
open $in, '<', "$top/RUNS/$RunName/$trait"."_Statistics.txt";
my $head = (<$in>);

print "Selecting Significant Results at FDR p-value < $bjh_alpha or average Bonferroni p-value < $avg_p (alpha=$bonf_alpha)!\n";
print "Also printing results where p < $p_annotate for completeness.\n";

while (<$in>) { #Prints Significant Hits to output file

	chomp $_;
	my @names = split("\t", $_);
	my $SNP=$names[0];
	my $Chrom=$names[1];
	my $BP=$names[2];
	my $pval=$names[3];
	my $fdrpval=$names[8];
	
	#if (($fdrpval < $bjh_alpha) || ($pval < $pcalc[$Chrom-1]) || ($pval < $p_annotate)) { #must be below BJH, bonferroni, AND user-specified threshold
	#if ($fdrpval < $bjh_alpha) { #must be below only BJH threshold
	if ($pval < $pcalc[$Chrom-1]) { #must be below only user threshold
	
		my $tmpos = join(",",$SNP,$Chrom,$BP,$pval,$fdrpval);
		chomp $tmpos;
		print $out "$tmpos\n";
		
	}
	
}

close $out;
close $in;

open $in, '<', "$top/RUNS/$RunName/$trait"."_SignificantHits";
open my $finalarabidopsis, '>', "$top/RUNS/$RunName/$trait"."_Arabidopsis.txt";
open my $finalrice, '>', "$top/RUNS/$RunName/$trait"."_Rice.txt";

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

print "Final Report Assembled!\n";

unless ($debug == 1) {

	unlink "$top/RUNS/$RunName/$trait"."_StatisticsTMP.txt";
	
}

mkdir "$top/RUNS/$RunName/INTERMEDIATES";
move("$top/RUNS/$RunName/$trait"."_Statistics.txt", "$top/RUNS/$RunName/INTERMEDIATES/$trait"."_Statistics.txt");
move("$top/RUNS/$RunName/$trait"."_Effects.txt", "$top/RUNS/$RunName/INTERMEDIATES/$trait"."_Effects.txt");
move("$top/RUNS/$RunName/$trait"."_SignificantHits", "$top/RUNS/$RunName/INTERMEDIATES/$trait"."_SignificantHits.txt");

open $out, '>', "$top/RUNS/$RunName/INTERMEDIATES/$trait"."_Significance_Thresholds.txt";
print $out "Bonferroni Alpha: $bonf_alpha\nFDR-corrected p-value Threshold: $bjh_alpha\nAverage Bonferroni Threshold p-value: $avg_p\nAverage Bonferroni Threshold Y-value: $avg_Y\n";
close $out;

print "Cleanup complete!\n";

}

print "GWAS PIPELINE COMPLETE!\n";
print "Press ENTER to Exit.\n>";
my $DONE = <>;