#script for producing generic plots of genome architecture from .gff3 file and .gff3 file-derived gene position list

#configuration options are:
#	.gff3 file for target genome: 40
#	feature coloring for the plot: starts line 99
#	types of featuers to plot: 152-154

#STRUCTURE
#SD3/PLOTTER/
#├── SD3-S4.pl #script to produce genome architecture plots for areas of DNA between given genes based on .gff3 file
#├── input *Input.txt file(s) containing beginning and end genes to include in plot, each on a separate line
#├── output *png.txt file(s) #plot
#└── COREFILES/
#	 ├── GenePositions.txt #list of gene positions (SD1-S6)
#	 └── Zea_mays.B73_RefGen_v4.45.gff3 #gene information for target genome

use strict;
use warnings;
use Cwd qw(abs_path);
use File::Copy;
use File::Path;
use File::Basename;

my $type = fileparse_set_fstype("Unix"); #Allows spaces in file paths
my $dirname = dirname(__FILE__);
my $top = defined($dirname) ? $dirname : '.';

my %chrom;
open my $mid, '<', "$top/COREFILES/GenePositions.txt";
while (<$mid>) {

	chomp $_;
	my @columns = split("\t", $_);
	$chrom{$columns[2]} = $columns[0];

}

close $mid;

my $GFF = "$top/COREFILES/Zea_mays.B73_RefGen_v4.45.gff3";
my @Input = glob("$top/*Input.txt");

open my $in, '<', $GFF;
$/ = '###';

my %start;
my %stop;
while (<$in>) {

	if ($_ =~ m/ID=gene:(.*?);/) {
	
		my $gene = $1;
		my @columns = split("\t", $_);
		
		$start{$gene} = $columns[3];
		$stop{$gene} = $columns[4];
	
	}

}

close $in;
$/ = "\n";

for my $file (@Input) {

$file =~ m/(.*?)_Input.txt/;
my $trait = $1;

open my $in, '<', $file;

my $first = <$in>;
chomp $first;
my $second = <$in>;
chomp $second;
close $in;

my @coordset = ($start{$first},$stop{$first},$start{$second},$stop{$second});
my @sorted_coords = sort { $a <=> $b } @coordset;

my $Start = shift(@sorted_coords);
my $Stop = pop(@sorted_coords);
my $Length = $Stop - $Start;
my $Chrom = $chrom{$first};

my $color = '"#ffd700"';
open my $py, '>', "python.py";
print $py "import numpy as np\nimport pandas as pd\nfrom pandas import DataFrame\nimport matplotlib.pyplot as plt\nfrom dna_features_viewer import (\n\tGraphicFeature,\n\tGraphicRecord,\n)\n\n";
print $py "font12 = {'family' : 'Arial',\n\t'color': 'black',\n\t'weight': 'bold',\n\t'size': 12,\n\t}\n\nfont14 = {'family' : 'Arial',\n\t'color': 'black',\n\t'weight': 'bold',\n\t'size': 14,\n\t}\n\n";
print $py "plt.rcParams[\"axes.linewidth\"] = 1.5\nfig, ax1 = plt.subplots(\n\t1, figsize = (12, 9)\n)\n\nfeatures = [\n";

my $sandwichtop = "\tGraphicFeature(\n\t\t";
my $sandwichbottom = "\n\t),\n";
my @annotations;
my %unique;

my %color = ( 

	#unusual features
	RNase_MRP_RNA => '"#FFFFFF"',
	SRP_RNA => '"#FFFFFF"',
	pseudogene => '"#FFFFFF"',
	ncRNA => '"#FFFFFF"',
	ncRNA_gene => '"#FFFFFF"',

	#RNA features
	rRNA => '"#FF0000"', 
	tRNA => '"#FF0000"',
	lnc_RNA => '"#FF0000"',
	pre_miRNA => '"#FF0000"',
	miRNA => '"#FF0000"',
	snRNA => '"#FF0000"',
	snoRNA => '"#FF0000"',

	#highly relevant features
	five_prime_UTR => '"#00FF00"',
	gene => '"#FFFFFF"',
	exon =>'"#0000FF"',
	CDS => '"#007FFF"', #not plotted
	three_prime_UTR =>  '"#7F00FF"',
	
	mRNA => '"#000000"', #not plotted

);

sub annotate {

	my @tannotations;
	my $entry = $_[0];
	my @lines = split("\n", $entry);
	my @real_lines;
	
	for my $line (@lines) { 
	
		unless ($line !~ m/^\d/) { push @real_lines, $line; }
	
	}
	
	for my $line (@real_lines) {
	
		my @columns = split("\t", $line);
		my $type = $columns[2];
		my $beginning = $columns[3];
		my $end = $columns[4];
		my $strand = $columns[6];
		my $annot = $columns[8];
		my $name;
		
		$annot =~ m/^.*?:(.*?);/;
		$name = $1;
		
		#if ($type =~ /gene|UTR|exon|CDS|[^m]RNA/) {
		#if ($type =~ /gene|UTR|exon|[^m]RNA/) {
		if ($type =~ /gene|exon|[^m]RNA/) {
		
			my $ID = "$beginning"."$end"."$type";
			
			unless(exists($unique{$ID})) {
		
				$unique{$ID} = 1;
				
				if ($type =~ /gene/) {
				
					my $label = "$type: "."$name";
					my $string = "$sandwichtop"."start=$beginning, end = $end, strand=$strand"."1, color=$color{$type}, label=\"$label\", linewidth = 1.5, fontdict = font12"."$sandwichbottom";
					push @tannotations, $string;
				
				} elsif ($type =~ /[^m]RNA/) {
				
					my $label = "$type: "."$name";
					my $string = "$sandwichtop"."start=$beginning, end = $end, strand=$strand"."1, color=$color{$type}, label=\"$label\", linewidth = 1.5, fontdict = font12"."$sandwichbottom";
					push @tannotations, $string;
					
				} else {
				
					my $string = "$sandwichtop"."start=$beginning, end = $end, strand=$strand"."1, color=$color{$type}, linewidth = 1.5"."$sandwichbottom";
					push @tannotations, $string;
					
				}
				
			}
		
		}
	
	}

	return @tannotations;

}

open $in, '<', $GFF;
$/ = '###';

my $flag = 0;
while (<$in>) {

	if ($_ =~ m/ID=gene:(.*?);/) {
	
		chomp $_;
		my $gene = $1;
		
		if ($flag == 1) {
		
			push @annotations, annotate($_);
		
		}
		
		if ($gene eq $second) {
		
			push @annotations, annotate($_);
			last;
		
		}
		
		
		if ($gene eq $first) {
			
			$flag = 1;
			push @annotations, annotate($_);
			
		}
	
	}

}

close $in;
$/ = "\n";

for my $entry (@annotations) {

	print $py $entry;

}

my $buffer_amt = $Length/100;
my $buffernear = $Start-$buffer_amt;
my $bufferfar = $Length+(2*$buffer_amt);

print $py "]\n\nrecord = GraphicRecord(sequence_length=$bufferfar,first_index=$buffernear, features=features)\nrecord.plot(ax=ax1, with_ruler=False, strand_in_label_threshold=4)\n\n";
print $py "ax1.set_xlabel(\"Chromosome $Chrom Position, bp\", fontdict=font14)\nax1.tick_params(width='2')\nax1.minorticks_off()\nplt.xticks(fontname='arial',fontweight='bold',fontsize=12)\n\n";
print $py "fig.savefig(\"$trait.pdf\", dpi=400, bbox_inches='tight')";
close $py;

system("python.py");
unlink "python.py";

}