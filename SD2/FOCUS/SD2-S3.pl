#configuration options start at line 52 for:
#	threshold options for selecting candidate significant markers. which of these options to use is defined at lines 617-618 and 680-682
#	method for assigning markers to genes-- specified LD window or n closest genes to marker

#more in-depth customization:
#	number of markers in datasets 2.7 and 3.2.1 at lines 65-66 and 328-337
#	GAPIT model defined at lines 555, together with surrounding script. For example, enable compression by setting the group.from parameter to 0, and/or by adding a group.to parameter to adjust the groupings.

#STRUCTURE
#SD2/FOCUS/
#├── SD2-S3.pl #script to execute COMPILE
#├── input *.txt file(s) containing phenotype data
#├── RUNS/
#│	 ├── Example.txt #contains example input data in correct format, for FOCUS also specify regions to execute GWAS in
#│	 └── RESULTS_FOLDERS #one per input data file 
#│			├── *_DATA.txt #copy of original input phenotype data
#│			└── RESULTS_SUBFOLDERS #one per chromosome region specified in input file
#│			 	 ├── INTERMEDIATES/ #raw GAPIT output files, including some not used by COMPILE (e.g. reports of marker effect size)
#│				 ├── output *_Significance_Thresholds.txt #report of significance thresholds used for candidate selection
#│				 ├── output *_Manhattan.pdf #Manhattan plots of results
#│				 ├── output *_QQ.pdf #quantile-quantile plot of results
#│				 ├── output *_Arabidopsis #candidate gene report matching maize to arabidopsis genes
#│				 └── output *_Rice #candidate gene report matching maize to rice genes
#└── COREFILES/
#	 ├── GAPIT/ #contains files used by GAPIT to execute GWAS
#	 └── GENOME/
#		 ├── 2.7/ #marker data for Goodman 282 low-density markers
#		 │ 	 ├── 1.GD.txt - 10.GD.txt, 1.GM.txt - 10.GM.txt #marker data in GAPIT numerical format (SD1-S2, SD1-S3)
#		 │ 	 ├── 1k.txt - 10k.txt #kinship matrix data (SD1-S4)
#		 │ 	 └── taxa.txt #list of taxa
#		 ├── 3.2.1/ #marker data for Goodman 282 high-density markers
#		 │ 	 ├── 1.GD.txt - 10.GD.txt, 1.GM.txt - 10.GM.txt #marker data in GAPIT numerical format (SD1-S2, SD1-S3)
#		 │ 	 ├── 1k.txt - 10k.txt #kinship matrix data (SD1-S4)
#		 │ 	 └── taxa.txt #list of taxa
#		 ├── GeneNameList.gff3 #filtered gene info list (SD1-S5)
#		 ├── GenePositions.txt #list of gene positions (SD1-S6)
#		 ├── KnownGenes.txt #list of known genes in maize
#		 ├── taxa.txt #list of taxa
#		 ├── Nearest.txt #atlas of nearest ten genes to each marker (SD1-S7)
#		 └── Rice.txt, Arabidopsis.txt #sequence similarity databases relating maize to rice and arabidopsis (SD1-S8, SD1-S9)

use strict;
use warnings;
use File::Basename;
use Cwd qw(abs_path);
use File::Copy;
use File::Path;
use List::Util qw(sum);
use POSIX qw(ceil);
use Number::Closest::NonOO qw(find_closest_number find_farthest_number);

my $debug = 0; #If 1, doesn't delete temporary files or directories
my $bjh_alpha = 0.1; #Significance threshold for BJH FDR p-val
my $bonf_alpha = 0.1; #Significance threshold for Bonferroni alpha
my $p_annotate = 0.0001; #Significance threshold to annotate markers anyway

my $use_ld = 0; #whether to use number of close genes or LD window. $number_of_close_genes option will still function, e.g. if >1 gene within the LD window, so set appropriately.
my $number_of_close_genes = 10; #number of genes near each marker to annotate (max 10)
my $ld_window = 10000; #window in bp around each marker to find genes to annotate

my $type = fileparse_set_fstype("Unix"); #Allows spaces in file paths
my $dirname = dirname(__FILE__);
my $top = defined($dirname) ? $dirname : '.';

my @num_markers_27 = (49384,37798,36146,29340,36591,24611,26431,26111,23166,20905);
my @num_markers_321 = (3165069,2498723,2436629,2374196,2014626,1609935,1603364,1679038,1702487,1466141);

my @bonferroni_27;
my @pcalc_27;

for my $nm (@num_markers_27) {

	push @pcalc_27, $bonf_alpha/$nm; #Y-values for significance thresholds on plot
	my $x = -(log($bonf_alpha/$nm)/log(10));
	push @bonferroni_27, $x;
	
}

my @bonferroni_321;
my @pcalc_321;

for my $nm (@num_markers_321) {

	push @pcalc_321, $bonf_alpha/$nm;
	my $x = -(log($bonf_alpha/$nm)/log(10));
	push @bonferroni_321, $x;
	
}

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
	$mappings{$columns[0]} = join("\t",@columns[1..10]);

}

close $map;

my %known;
open my $knownfile, '<', "$top/COREFILES/Genome/KnownGenes.txt";
my $header = <$knownfile>;
while (<$knownfile>) {

	chomp $_;
	my @columns = split("\t", $_);
	$known{$columns[0]} = join('/',@columns[1..2]);

}

close $knownfile;

my %rice;
open my $ricefile, '<', "$top/COREFILES/Genome/Rice.txt";
my $header1 = <$ricefile>;
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

my @AoH;
open my $translate, '<', "$top/COREFILES/Genome/GenePositions.txt";

while (<$translate>) {

	chomp $_;
	my @columns = split("\t", $_);
	$AoH[$columns[0]]{$columns[1]} = $columns[2];
	
}

close $translate;

open my $RPath, '<', "$top/COREFILES/R_Installation.txt";
my $R = <$RPath>;
close $RPath;

print "Please place copies of your input data in .txt format in the base Pipeline folder.\nThe Pipeline will run them sequentially.\nData files will be renamed and moved to the run folder in /RUNS/.\n";
print "Please ensure the input file is in the format in /RUNS/Template.txt. Coordinates are accepted in units of megabases.\nPress any key to continue.\n>";
my $GO = <>;

my @datafiles = glob("$top/*.txt");
for my $file (@datafiles) { #Executes the Pipeline for every set of trait data in the folder

	opendir my $dh, "$top/RUNS"; #Gets the name of every existing directory in RUNS
	my @existing_runs = grep {-d "$top/RUNS/$_" && ! /^\.{1,2}$/} readdir($dh);
	closedir $dh;

	open my $traitfile, '<', "$file"; #Gets the first line of the input file
	my $firstline = <$traitfile>;

	my @queries;
	my @phenotypes;

	my $trait;
	while (<$traitfile>) {

		chomp $_;
		if ($_  =~ m/Taxa\t(.*?)$/) {
	
			$trait = $1;
			while (<$traitfile>) {
		
				push(@phenotypes, $_);
			
			}
		
		}
	
	else { push(@queries, $_); }
	
	}

close $traitfile;
chomp $trait;
print "Working on queries for $trait!\n";
chomp @queries;
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
my $RunDataName = join('_', $trait,'DATA');
copy("$file", "$top/RUNS/$RunName/$RunDataName.txt");

my @missing;
chomp @phenotypes;
open my $phenos, '>', "$top/RUNS/$RunName/$RunDataName.txt";
print $phenos "Taxa\t$trait\n";
foreach(@phenotypes) { 

	print $phenos "$_\n";
	if ($_ =~ m/NaN/) { push @missing, (split("\t", $_))[0]; }
	
}

close $phenos;
chomp @missing;

my @versions = ("2.7","3.2.1");

for my $version (@versions) {
print "\tWorking on Version $version...\n";

open my $reference, '<', "$top/COREFILES/GENOME/$version/taxa.txt"; #Nab the taxa list
my @columns = split("\t", <$reference>); #Get an array of taxa names in the marker data file
close $reference;
chomp @columns;

my @keep; #Defines an array of taxa names to keep

for my $j (0..$#columns) { #Pushes only taxa for which there is data to @keep by looping through @missing for each value in @columns
	
	my $flg = 0;
	for my $k (0..$#missing) {
		
		if ($missing[$k] eq $columns[$j]) { $flg = 1; } #If the taxon is missing, set the flag to 0 so it doesn't get pushed to @keep
	
	}
	
	if ($flg == 0) { push @keep, $columns[$j]; }

}

for my $q (@queries) {

my $input_type;
my $position1;
my $position2;
my @local_keep = @keep;

my $bon;
my $threshold;
my @final_markers;
my $compression_taxa;

my $chromosome = (split("\t", $q))[0];
	next unless ($chromosome =~ /[\d]/);

if ($version eq '2.7') { 

	$bon = $pcalc_27[$chromosome-1];
	$threshold = $bonferroni_27[$chromosome-1];
	@final_markers = (0,306971061,244417267,235520195,246943212,223706202,173351536,182129497,181044202,159696575,150890390); 
	$compression_taxa = 282;

}

elsif ($version eq '3.2.1') { 

	$bon = $pcalc_321[$chromosome-1];
	$threshold = $bonferroni_321[$chromosome-1];
	@final_markers = (0,306970776,244440046,235653210,246967124,223706459,173377089,181718895,181046086,159687194,150930636); 
	$compression_taxa = 279;

}

my $posinput = (split("\t", $q))[1];
if ($posinput !~ /-/) { $input_type = "Single" }
else { $input_type = "Range" }

if ($input_type eq "Single") {

	my $suffix_type;
	my $position;

	$position = $posinput*1000000;

	if ($position >= $final_markers[$chromosome]) {
	
		$position1 = $final_markers[$chromosome]-10000001;
		$position2 = $final_markers[$chromosome]-1; 
	
	}

	else { 
	
		$position1 = $position-5000000;
		$position2 = $position+5000000;
	
	}
	
}

if ($input_type eq "Range") {

	my $pos1 = (split("-", $posinput))[0];
	chomp $pos1;
	$position1 = $pos1*1000000;

	if ($position1 >= $final_markers[$chromosome]) {
	
	$position1 = $final_markers[$chromosome]-1; 
	
	}

	my $pos2 = (split("-", $posinput))[1];
	chomp $pos2;
	$position2 = $pos2*1000000;

	if ($position2 >= $final_markers[$chromosome]) {
	
	$position2 = $final_markers[$chromosome]-1; 
	
	}

	if ($position1 == $position2) { $position1 = $position1-10000000; }

	my @positions = ($position1,$position2);
	my @sorted_positions = sort { $a <=> $b } @positions;
	$position1 = $sorted_positions[0];
	$position2 = $sorted_positions[1];

}

mkdir "$top/RUNS/$RunName/$version-$chromosome-$posinput", 0755;
mkdir "$top/RUNS/$RunName/$version-$chromosome-$posinput/TMP", 0755;

my $pos1MB = $position1/1000000;
my $pos2MB = $position2/1000000;

open my $manhattan, '<', "$top/COREFILES/GAPIT/manhattan.txt";
my $file_content = do { local $/; <$manhattan> };
open my $local_manhattan, '>', "$top/RUNS/$RunName/$version-$chromosome-$posinput/TMP/manhattan.txt";
$file_content =~ s/insertxmin/xmin = $pos1MB/;
$file_content =~ s/insertxmax/xmax = $pos2MB/;
print $local_manhattan $file_content;
close $manhattan;
close $local_manhattan;

open my $index, '<', "COREFILES/GENOME/$version/$chromosome.GM.txt";
my $header = <$index>;

my $index_first;
my $index_last;
my $counter = 0;

while (<$index>) {

	my $position = (split("\t", $_))[2];
	chomp $position;

	if ($position1 >= $position) { 
	
		$counter++;
		next; 
		
	}
	
	elsif ($position1 <= $position) {	
	
		$index_first = $counter;
		$counter++;
		last;
		
	}
	
}

while (<$index>) {

	my $position = (split("\t", $_))[2];
	chomp $position;

	if ($position2 >= $position) { 
	
		$counter++;
		next; 
		
	}
	
	elsif ($position2 <= $position) {	
	
		$index_last = $counter;
		$counter++;
		last;
		
	}
	
}

close $index;

my $taxa = scalar(@keep); #Number of taxa is equal to the number of kept columns
my %to_keep = map { $_ => 1 } @keep;

print "\t\tExecuting MLM via GAPIT for Chromosome $chromosome & $posinput ($taxa Individuals)...\n";

open my $in, '<', "COREFILES/GENOME/$version/$chromosome.GD.txt";
open my $out, '>', "$top/RUNS/$RunName/$version-$chromosome-$posinput/TMP/$chromosome.GD.txt";

my $firstline = <$in>;
my @firstlinecolumns = split("\t", $firstline);
my @headeroutput;

push @headeroutput, $firstlinecolumns[0];

for my $headernum ($index_first+1..$index_last+1) {
		
	push @headeroutput, $firstlinecolumns[$headernum];
		
}

my $GDHead = join("\t", @headeroutput); 
print $out "$GDHead\n"; 

while (<$in>) {

	chomp $_;
	
	
	my $taxon = (split("\t", $_))[0];
	if(exists($to_keep{$taxon})) { 
	
		my @cols = split("\t", $_);
		
		my @output;
		push @output, $cols[0];
		for my $number ($index_first+1..$index_last+1) {
		
			push @output, $cols[$number];
		
		}
		
		my $outline = join("\t", @output); 
		print $out "$outline\n"; 
	
	}	
	
}

open $in, '<', "COREFILES/GENOME/$version/$chromosome.GM.txt";
open $out, '>', "$top/RUNS/$RunName/$version-$chromosome-$posinput/TMP/$chromosome.GM.txt";
$header = <$in>;
print $out $header;

my @GM_lines = <$in>;
close $in;

for my $number ($index_first..$index_last) {

	print $out $GM_lines[$number];

}

close $out;

my $k = "$chromosome"."k";
my $I = $chromosome;

open my $rfile, '>', "$top/RUNS/$RunName/$version-$chromosome-$posinput/TMP/R1.txt";
my $core = dirname(abs_path($0));

$core =~ s|\\|/|g;
$core =~ s|(.*?)/Pipeline.pl|$1|g;

print $rfile qq`library(compiler)
source("$core/COREFILES/GAPIT/gapit_functions.txt")
source("$core/COREFILES/GAPIT/emma.txt")
source("$core/COREFILES/GAPIT/mt.txt")
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(scatterplot3d)
myGM <- read.table("$core/RUNS/$RunName/$version-$chromosome-$posinput/TMP/$chromosome.GM.txt", head=TRUE)
myGD <- read.table("$core/RUNS/$RunName/$version-$chromosome-$posinput/TMP/$chromosome.GD.txt", head=TRUE)
myKI <- read.table("$core/COREFILES/GENOME/$version/$k.txt", header=FALSE)
myY <- read.table("$core/RUNS/$RunName/$RunDataName.txt", head=TRUE)
setwd("$core/RUNS/$RunName/$version-$chromosome-$posinput/TMP")
myGAPIT <- GAPIT(Y=myY, GD=myGD, GM=myGM, KI=myKI, group.from=$compression_taxa, group.to=$compression_taxa, Geno.View.output=FALSE)
q()
`;

close $rfile;
my $command = "$R/Rscript $top/RUNS/$RunName/$version-$chromosome-$posinput/TMP/R1.txt > nul 2>\&1";
system("$command");

#Compiles GWAS Results
open $rfile, '>', "$top/RUNS/$RunName/$version-$chromosome-$posinput/TMP/R2.txt";
print $rfile qq`library(compiler)
source("$core/COREFILES/GAPIT/gapit_functions.txt")
source("$core/COREFILES/GAPIT/mt.txt")
source("$core/RUNS/$RunName/$version-$chromosome-$posinput/TMP/manhattan.txt")
library(gplots)
library(genetics)
library(EMMREML)
library("scatterplot3d")
library(MASS)
library(manhattanly)
library(plotly)
name <- "$trait"
GWAS.Results <- read.csv("$core/RUNS/$RunName/$version-$chromosome-$posinput/TMP/GAPIT..$trait.GWAS.Results.csv",head=TRUE)
Effect.Estimates <- read.csv("$core/RUNS/$RunName/$version-$chromosome-$posinput/TMP/GAPIT..$trait.Allelic_Effect_Estimates.csv", head=TRUE)
GWAS.Results <- rbind(GWAS.Results)
Effect.Estimates <- rbind(Effect.Estimates)
GWAS.Results <- GWAS.Results[,-ncol(GWAS.Results)]
setwd("$core/RUNS/$RunName/$version-$chromosome-$posinput/TMP")
Conduct.FDR <- GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(PWI=GWAS.Results, FDR.Rate=$bjh_alpha, FDR.Procedure="BH")
GWAS.Results.FDR <- Conduct.FDR\$PWIP
setwd("$core/RUNS/$RunName/$version-$chromosome-$posinput/TMP")
write.table(GWAS.Results.FDR, paste("GAPIT.", name, ".GWAS.Results.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(Effect.Estimates, paste("GAPIT.", name, ".Allelic_Effect_Estimates.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
q()
`;

close $rfile;
$command = "$R/Rscript $top/RUNS/$RunName/$version-$chromosome-$posinput/TMP/R2.txt > nul 2>\&1";
system("$command");

move("$top/RUNS/$RunName/$version-$chromosome-$posinput/TMP/GAPIT.$trait.Allelic_Effect_Estimates.txt", "$top/RUNS/$RunName/$version-$chromosome-$posinput/$trait"."_Effects.txt");
move("$top/RUNS/$RunName/$version-$chromosome-$posinput/TMP/GAPIT.$trait.GWAS.Results.txt", "$top/RUNS/$RunName/$version-$chromosome-$posinput/$trait"."_StatisticsTMP.txt");

open $in, '<', "$top/RUNS/$RunName/$version-$chromosome-$posinput/$trait"."_StatisticsTMP.txt"; #Creates a new statistics file sorted by corrected FDR P-value
my $statsheader = <$in>;
my @lines = <$in>;
chomp @lines;
close $in;

my @ordered_lines = sort { (split("\t", $a))[1] <=> (split("\t", $b))[1] || (split("\t", $a))[2] <=> (split("\t", $b))[2] } @lines;

open $out, '>', "$top/RUNS/$RunName/$version-$chromosome-$posinput/$trait"."_Statistics.txt";
print $out $statsheader;
for (@ordered_lines) { print $out "$_\n"; }
close $out;

my @BJH; #Pushes to @BJH FDR P-values below $alpha, these hits are colored in the final Manhattan plot
open $in,  '<',  "$top/RUNS/$RunName/$version-$chromosome-$posinput/$trait"."_Statistics.txt";
my $discard = (<$in>);

while (<$in>) {

	#if ((split("\t", $_))[8] < $bjh_alpha) { #bjh-corrected p-value below BJH threshold
	if ((split("\t", $_))[3] < $bon) { #uncorrected p-value below Bonferroni threshold
	
		push @BJH, (split("\t", $_))[0]; 
		
	}
	
}

close $in;
my $R_HITS = join('","', @BJH);

#Assimilate the final interactive Manhattan plots

open $rfile, '>', "$top/RUNS/$RunName/$version-$chromosome-$posinput/TMP/R3.txt";
print $rfile qq`library(compiler)
setwd("$top/RUNS/$RunName/$version-$chromosome-$posinput")
source("$core/COREFILES/GAPIT/gapit_functions.txt")
source("$core/COREFILES/GAPIT/mt.txt")
library(gplots)
library(genetics)
library(EMMREML)
library("scatterplot3d")
library(MASS)
source("$core/RUNS/$RunName/$version-$chromosome-$posinput/TMP/manhattan.txt")
library(manhattanly)
library(plotly)
name <- "$trait"
GWAS.Results <- read.table(paste("$core/RUNS/$RunName/$version-$chromosome-$posinput/$trait", "_StatisticsTMP.txt", sep = ""), head=TRUE)
SNPsOfInterest <- c("$R_HITS")
GAPIT.QQ(P.values = GWAS.Results[,4], name.of.trait = name,DPP=50000)
manhattan(GWAS.Results, genomewideline=$threshold, chr="Chromosome", bp="Position", p="P.value", snp="SNP", col=c("black", "gold3"), main = "Manhattan Plot for $trait", highlight = SNPsOfInterest)
GR <- manhattanr(GWAS.Results, chr = "Chromosome", bp = "Position", p = "P.value")
htmlwidgets::saveWidget(manhattanly(subset(GR[["data"]], CHR %in% $chromosome), genomewideline=$threshold, title = "Manhattan Plot for $trait", col=c("black")), "Chr_$chromosome.html")
q()
`;

close $rfile;
$command = "$R/Rscript $top/RUNS/$RunName/$version-$chromosome-$posinput/TMP/R3.txt > nul 2>\&1";
system("$command");

#Finish sorting out GAPIT outputs
move("$top/RUNS/$RunName/$version-$chromosome-$posinput/Rplots.pdf", "$top/RUNS/$RunName/$version-$chromosome-$posinput/$trait"."_Manhattan.pdf");
move("$top/RUNS/$RunName/$version-$chromosome-$posinput/GAPIT.$trait.QQ-Plot.pdf", "$top/RUNS/$RunName/$version-$chromosome-$posinput/$trait"."_QQ.pdf");

unless ($debug == 1) { rmtree("$top/RUNS/$RunName/$version-$chromosome-$posinput/TMP"); }

my @Hi_Mom;

open $out, '>', "$top/RUNS/$RunName/$version-$chromosome-$posinput/$trait"."_SignificantHits";
open $in, '<', "$top/RUNS/$RunName/$version-$chromosome-$posinput/$trait"."_Statistics.txt";
my $head = (<$in>);

while (<$in>) { #Prints Significant Hits to output file

	chomp $_;
	my @names = split("\t", $_);
	my $SNP=$names[0];
	my $Chrom=$names[1];
	my $BP=$names[2];
	my $pval=$names[3];
	my $fdrpval=$names[8];
	
	#if (($fdrpval < $bjh_alpha) || ($pval < $bon) || ($pval < $p_annotate)) { #must be below BJH, bonferroni, AND user-specified threshold
	#if ($fdrpval < $bjh_alpha) { #must be below only BJH threshold
	if ($pval < $bon) { #must be below only user threshold
	
		my $tmpos = join(',',$SNP,$Chrom,$BP,$pval,$fdrpval);
		chomp $tmpos;
		print $out "$tmpos\n";
		
	}
	
}

close $out;
close $in;


#Find maize gene names from significant hit marker positions
open $in, '<', "$top/RUNS/$RunName/$version-$chromosome-$posinput/$trait"."_SignificantHits";
open my $finalarabidopsis, '>', "$top/RUNS/$RunName/$version-$chromosome-$posinput/$trait"."_Arabidopsis.txt";
open my $finalrice, '>', "$top/RUNS/$RunName/$version-$chromosome-$posinput/$trait"."_Rice.txt";

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
	
	my @genelist;

	if ($version eq "2.7") { @genelist = split("\t", $mappings{$SNP}); }
	elsif ($version eq "3.2.1") {
	
		my @coords = keys(%{$AoH[$chromosome]});
		@coords = sort { $a <=> $b } @coords;
	
		my $closest = find_closest_number(number=>$posit, numbers=>\@coords, items => 10);
		my @values = @{$closest};
		@values = sort @values;

		for my $k (0..9) {
		
			push @genelist, $AoH[$chromosome]{$values[$k]};
	
		}
	
	
	}
	
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

unless ($debug == 1) { unlink "$top/RUNS/$RunName/$version-$chromosome-$posinput/$trait"."_StatisticsTMP.txt"; }

mkdir "$top/RUNS/$RunName/$version-$chromosome-$posinput/INTERMEDIATES", 0755;
move("$top/RUNS/$RunName/$version-$chromosome-$posinput/$trait"."_SignificantHits", "$top/RUNS/$RunName/$version-$chromosome-$posinput/INTERMEDIATES/$trait"."_SignificantHits.txt");
move("$top/RUNS/$RunName/$version-$chromosome-$posinput/$trait"."_Statistics.txt", "$top/RUNS/$RunName/$version-$chromosome-$posinput/INTERMEDIATES/$trait"."_Statistics.txt");
move("$top/RUNS/$RunName/$version-$chromosome-$posinput/$trait"."_Effects.txt", "$top/RUNS/$RunName/$version-$chromosome-$posinput/INTERMEDIATES/$trait"."_Effects.txt");

open $out, '>', "$top/Runs/$RunName/$version-$chromosome-$posinput/$trait"."_Significance_Thresholds.txt";
print $out "Bonferroni Alpha: $bonf_alpha\nFDR-corrected p-value Threshold: $bjh_alpha\nAverage Bonferroni Threshold p-value: $bon\nAverage Bonferroni Threshold Y-value: $bon\n";
close $out;

}

}

}

print "GWAS PIPELINE COMPLETE!\n";
print "Press ENTER to Exit.\n>";
my $DONE = <>;
