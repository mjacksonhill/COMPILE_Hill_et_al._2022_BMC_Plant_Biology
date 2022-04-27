#identify the ten nearest genes (from marker .gm.txt files in GAPIT numerical format to each specified gene in GenePositions.txtmarker set
#requires the marker GM files used by GAPIT (product of SD1-S3.pl) and GenePositions.txt (product of SD1-S6.pl)

#number of nearest genes specified on lines 48, 52
#number of chromosomes specified on line 32

#STRUCTURE
#SD1/
#├── SD1-S7.pl #script to run after SD1-S6.pl to generate list of ten nearest genes to each marker
#├── Input 1.GM.txt - ##.GM.txt #genotype files for panel in GAPIT numerical format (e.g. from SD1-S3.pl)
#├── GenePositions.txt file from SD1-S6.pl
#└── Nearest.txt output file

use strict;
use warnings;
use Number::Closest;

my @AoH;
open my $translate, '<', 'GenePositions.txt';

while (<$translate>) {

	chomp $_;
	my @columns = split("\t", $_);
	$AoH[$columns[0]]{$columns[1]} = $columns[2];
	
}

close $translate;
open my $output, '>', "Nearest.txt";

for my $j (1..10) {

	open my $markers, '<', "$j.GM.txt";
	my $header = <$markers>;
	
	my @coords = keys(%{$AoH[$j]});
	@coords = sort { $a <=> $b } @coords;
	
	while (<$markers>) {
		
		chomp $_;
		my @columns = split("\t", $_);
		my $name = $columns[0];
		my $position = $columns[2];

		my $finder = Number::Closest->new(number => $position, numbers => \@coords) ;
		my @values = @{$finder->find(10)};
		@values = sort @values;
		
		my @temp;
		for my $k (0..9) {
		
			push @temp, $AoH[$j]{$values[$k]};
	
		}
		
		my $out = join("\t", @temp);
		my $outline = "$name\t$out";
		print $output "$outline\n";
	
	}
	
	close $markers;

	
}

close $output;