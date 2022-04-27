#convert the GeneNameList.gff3 file used by COMPILE (from SD1-S5.pl) into a list of gene positions for further use with SD1-S7.pl

#maize gene name format specified in regular expression on line 30

#STRUCTURE
#SD1/
#├── SD1-S6.pl #script to run to generate gene position list
#├── GeneNameList.gff3 input file (product of SD1-S5.pl)
#└── GenePositions.txt output file

use strict;
use warnings;

open my $in, '<', 'GeneNameList.gff3';
open my $out, '>', 'GenePositions.txt';

my %HoA;
my %Chromosomes;

while (<$in>) {

	chomp $_;
	my @columns = split("\t", $_);
	
	my $chromosome = $columns[0];
	my $start = $columns[2];
	my $stop = $columns[3];
	my $name = $columns[4];
	
	$name =~ s/(Zm.*?)_.{4}/$1/;
	push @{ $HoA{$name} }, $start, $stop;
	$Chromosomes{$name} = $chromosome;
	
}

close $in;
my @genes = keys %HoA;
my @outputs;

for my $gene (@genes) {
	
	my @positions = @{$HoA{$gene}};
	@positions = sort { $a <=> $b } @positions;

	my $chromosome = $Chromosomes{$gene};
	my $first = shift @positions;
	my $last = pop @positions;
	my $average = ($first+$last)/2;
	
	push @outputs, "$chromosome\t$average\t$gene\n";
	
}

my @sorted = sort {
    my ($aa, $bb) = map { (split /\t/)[0] } $a, $b;
	my ($cc, $dd) = map { (split /\t/)[1] } $a, $b;
    $aa <=> $bb
		or
	$cc <=> $dd;
} @outputs;


for my $line (@sorted) {

	print $out $line;

}

close $out;