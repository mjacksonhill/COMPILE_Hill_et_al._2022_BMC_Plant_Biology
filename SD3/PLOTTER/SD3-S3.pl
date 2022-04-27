#script for overlaying FOCUS dual Manhattan plots with maize local genome architecture

#configuration options are:
#	annotation colors: 32-33
#	.gff3 file for target genome: 51
#	feature coloring for the plot: starts line 200
#	types of featuers to plot: 253-255

#STRUCTURE
#SD3/PLOTTER/
#├── SD3-S3.pl #script to re-analyze GAPIT statistics files, generating new lists of significant markers
#├── input *_Input.txt file(s) containing beginning and end genes to include in plot, each on a separate line
#├── input *_Statistics_2.7.txt file(s) containing GAPIT results from low-density marker data
#├── input *_Statistics_3.2.1.txt file(s) containing GAPIT results from high-density marker data
#├── input *_SignificantHits_2.7.txt file(s) containing COMPILE-annotated significant genes from analysis of low-density marker data
#├── input *_SignificantHits_3.2.1.txt file(s) containing COMPILE-annotated significant genes from analysis of high-density marker data
#├── output *_png.txt file(s) #plot
#└── COREFILES/
#	 ├── GENOME/2.7/1.GM.txt - 10.GM.txt #Goodman panel marker data in GAPIT numerical format (SD1-S2, SD1-S3)
#	 ├── GENOME/3.2.1/1.GM.txt - 10.GM.txt #Goodman high-density panel marker data in GAPIT numerical format (SD1-S2, SD1-S3)
#	 ├── GenePositions.txt #list of gene positions (SD1-S6)
#	 └── Zea_mays.B73_RefGen_v4.45.gff3 #gene information for target genome

use strict;
use warnings;
use Cwd qw(abs_path);
use File::Copy;
use File::Path;
use File::Basename;
use POSIX;

my $basecolor = "black";
my $accentcolor = "red";

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
my @Input = glob("$top/*_Input.txt");

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

my %bonferroni27;
my %bonferroni321;
my $bonf_alpha= 0.1;

for my $version ('2.7', '3.2.1') {
	
	my %bonferroni;
	my @num_markers;
	for my $i (1..10) { #Get numbers of markers on each chromosome

		open my $file, '<', "$top/COREFILES/GENOME/$version/$i.GM.txt";
		while (<$file>) {}
		my $line = $.-1; #0-index takes care of -1 for empty line at end; another -1 for header
		push @num_markers, $line;
		close $file;

	}

	my $i = 1;
	
	for my $nm (@num_markers) {

		my $x = -(log($bonf_alpha/$nm)/log(10));
		$bonferroni{$i} = $x;
		$i++;
	
	}
	
	if ($version eq '2.7') { %bonferroni27 = %bonferroni; }
	else { %bonferroni321 = %bonferroni; }

}

for my $file (@Input) {

$file =~ m/\.\/(.*?)_Input.txt/;
my $trait = $1;

open my $in,  '<', $file;

my @one;
my @two;
$/ = "\n\n";
while (<$in>) {
	
	chomp $_;
	my $first = (split("\n", $_))[0];
	$first =~ s/\R//g;
	push @one, $first;
	my $second = (split("\n", $_))[1];
	$second =~ s/\R//g;
	push @two, $second;
	
}

close $in;
$/ = "\n";

for my $version ('2.7', '3.2.1') {
	
	$/ = undef;
	my $stats_string = "SNP	Chromosome	Position	P.value	maf	nobs	Rsquare.of.Model.without.SNP	Rsquare.of.Model.with.SNP	FDR_Adjusted_P-values\n";
	my $stats_replacement = "SNP\tChromosome\tposition\tpvalue\tx\ttaxa\tx\ty\tFDR.p.value\n";

	open $in,  '<', "$top/$trait"."_Statistics_$version.txt";
	my $file = <$in>;
	close $in;
	
	$file =~ s/$stats_string/$stats_replacement/;
	
	open my $out, '>', "$top/$trait"."_Statistics_$version.txt";
	print $out $file;
	close $out;
	
	my $sig_string = "SNP,Chromosome,position,pvalue,FDR.p.value\n";
	
	open $in,  '<', "$top/$trait"."_SignificantHits_$version.txt";
	$file = <$in>;
	close $in;
	
	unless ($file =~ m/$sig_string/) {
	
		open my $out, '>', "$top/$trait"."_SignificantHits_$version.txt";
		print $out $sig_string;
		print $out $file;
		close $out;
		
	}

}

for my $i (0..$#one) {
	
my $first = $one[$i];
my $second = $two[$i];

my @coordset = ($start{$first},$stop{$first},$start{$second},$stop{$second});
my @sorted_coords = sort { $a <=> $b } @coordset;

my $Start = shift(@sorted_coords);
my $Stop = pop(@sorted_coords);
my $Length = $Stop - $Start;
my $Chrom = $chrom{$first};

open $in, '<', $GFF;
$/ = '###';

my $color = '"#ffd700"';

my $pythonname = "python_$trait"."_$i.py";
open my $py, '>', $pythonname;
print $py "import numpy as np\nimport pandas as pd\nfrom pandas import DataFrame\nimport matplotlib.pyplot as plt\nimport matplotlib.ticker as tkr\nfrom dna_features_viewer import (\n\tGraphicFeature,\n\tGraphicRecord,\n)\n\n";
print $py "def numfmt(x,pos):\n\ts = '{}'.format(x / 1000000)\n\treturn s\n\nxfmt = tkr.FuncFormatter(numfmt)\n\n";
print $py "font12 = {'family' : 'Arial',\n\t'color': 'black',\n\t'weight': 'bold',\n\t'size': 12,\n\t}\n\nfont14 = {'family' : 'Arial',\n\t'color': 'black',\n\t'weight': 'bold',\n\t'size': 14,\n\t}\n\n";
print $py "plt.rcParams[\"axes.linewidth\"] = 1.5\nfig, (ax1, ax2, ax3) = plt.subplots(\n\t3, figsize = (12, 9), sharex=True, gridspec_kw={\"height_ratios\": [3, 1, 1], 'hspace': 0}\n)\n\nfeatures = [\n";

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

	my %localunique = %unique;
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
			
			unless((exists($unique{$ID}) || (exists($localunique{$ID})))) {
		
				$localunique{$ID} = 1;
				
				if ($type =~ /gene/) {
				
					unless($type =~ /ncRNA_gene/) {
						
						my $label = "$type: "."$name";
						my $string = "$sandwichtop"."start=$beginning, end = $end, strand=$strand"."1, color=$color{$type}, label=\"$label\", linewidth = 1.5, fontdict = font12"."$sandwichbottom";
						push @tannotations, $string;
					
					}
				
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
	
	return \@tannotations, \%localunique;

}

my $flag = 0;
while (<$in>) {

	if ($_ =~ m/ID=gene:(.*?);/) {
	
		chomp $_;
		my $gene = $1;
		
		if ($gene eq $second) {
		
			my ($array, $hash) = annotate($_);
			push @annotations, @$array;
			%unique = (%unique, %$hash);
			last;
		
		}
		
		if ($flag == 1) {
			
			my ($array, $hash) = annotate($_);
			push @annotations, @$array;
			%unique = (%unique, %$hash);
		
		}
		
		if ($gene eq $first) {
		
			$flag = 1;
			my ($array, $hash) = annotate($_);
			push @annotations, @$array;
			%unique = (%unique, %$hash);
			
		}
	
	}

}

close $in;

for my $entry (@annotations) {

	print $py $entry;

}

my $buffer_amt = $Length/100;
my $buffernear = $Start-$buffer_amt;
my $bufferfar = $Length+(2*$buffer_amt);

sub magic_number {
   my $n = shift;
   return (log($n)/log(10))-floor(log($n)/log(10));
}

my $number = magic_number($Length);

my $inc = 1;
if ($number >= magic_number(500)) { $inc = 10 }
elsif ($number >= magic_number(250)) { $inc = 5 }
elsif ($number >= magic_number(125)) { $inc = 2.5 }
my $multiple = (10 ** ((floor(log($Length)/log(10))) -1 )) * $inc;

print $py "]\n\nrecord = GraphicRecord(sequence_length=$bufferfar,first_index=$buffernear, features=features)\nrecord.plot(ax=ax1, with_ruler=False, strand_in_label_threshold=4)\n\n";
print $py "df = pd.read_csv(\"$trait"."_Statistics_2.7.txt\",sep='\\t',skiprows=(0),usecols=[2,3],header=(0))\ndf['minuslog10pval'] = -np.log10(df.pvalue)\ndf.plot.scatter(x='position', y='minuslog10pval', color='$basecolor', ax=ax2)\n\n";
print $py "df2 = pd.read_csv(\"$trait"."_SignificantHits_2.7.txt\",sep=',',skiprows=(0),usecols=[2,3],header=(0))\ndf2['minuslog10pval'] = -np.log10(df2.pvalue)\ndf2.plot.scatter(x='position', y='minuslog10pval', color='$accentcolor', ax=ax2)\n\n";
print $py "df3 = pd.read_csv(\"$trait"."_Statistics_3.2.1.txt\",sep='\\t',skiprows=(0),usecols=[2,3],header=(0))\ndf3['minuslog10pval'] = -np.log10(df3.pvalue)\ndf3.plot.scatter(x='position', y='minuslog10pval', color='$basecolor', ax=ax3)\n\n";
print $py "df4 = pd.read_csv(\"$trait"."_SignificantHits_3.2.1.txt\",sep=',',skiprows=(0),usecols=[2,3],header=(0))\ndf4['minuslog10pval'] = -np.log10(df4.pvalue)\ndf4.plot.scatter(x='position', y='minuslog10pval', color='$accentcolor', ax=ax3)\n\n";
print $py "ax3.xaxis.set_major_formatter(xfmt)\nax3.xaxis.set_major_locator(plt.MultipleLocator($multiple))\n\nax2.set_ylim(bottom=0)\nax3.set_ylim(bottom=0)\nax2.set_ylabel(\"-log10(p)\", fontdict=font14)\nax3.set_ylabel(\"-log10(p)\", fontdict=font14)\nax3.set_xlabel(\"Chromosome $Chrom Position, Mbp\", fontdict=font14)\nax2.axhline($bonferroni27{$Chrom}, color=\"$accentcolor\")\nax3.axhline($bonferroni321{$Chrom}, color=\"$accentcolor\")\n\n";
print $py "ax2.tick_params(width='2')\nax3.tick_params(width='2')\nax3.minorticks_off()\n\nplt.xticks(fontname='arial',fontweight='bold',fontsize=12)\nplt.yticks(fontname='arial',fontweight='bold',fontsize=12)\n\n";
print $py "axes5 = ax2.twinx()\naxes5.set_ylabel('Goodman 2.7', fontdict=font14,rotation=270, labelpad=15)\naxes5.yaxis.set_ticks([])\n\naxes6 = ax3.twinx()\naxes6.set_ylabel('Goodman 3.2.1', fontdict=font14,rotation=270, labelpad=15)\naxes6.yaxis.set_ticks([])\n\n";
print $py "plt.sca(ax2)\nplt.yticks(fontname='arial',fontweight='bold',fontsize=12)\nfig.savefig(\"$trait"."_$i.pdf\", dpi=400, bbox_inches='tight')";
close $py;

system("$pythonname");
unlink $pythonname;

}

}