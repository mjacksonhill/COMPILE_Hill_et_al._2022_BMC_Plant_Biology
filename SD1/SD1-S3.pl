#convert HapMap files to GAPIT numerical format files-- does not precisely follow rules as to reference alleles, but this is irrelevant for GWAS purposes.

#chromosome number specified on line 15
#missing data is imputed as the major allele at line 50; change as desired. GAPIT requires all data to be imputed, so options are major allele, heterozygote, or minor allele

#STRUCTURE
#SD1/
#├── SD1-S3.pl #script to perform conversion
#├── Output .GM.txt and .GD.txt file(s)
#└── Input 1.hmp.txt - ##.hmp.txt genotype files

use warnings;
use strict;

for my $i (1..10) {

	open my $infile, '<', "$i.hmp.txt";
	
	my %numerical;
	open my $GM, '>', "$i.GM.txt";
	print $GM "Name\tChromosome\tPosition\n";
	open my $GD, '>', "$i.GD.txt";
	
	my $header = <$infile>;
	my @taxa = split("\t", $header);
	@taxa = @taxa[11..$#taxa];
	chomp @taxa;
	
	my @markerids;
	
	while (<$infile>) {
	
		my @data = split("\t", $_);
		chomp @data;
		
		my $marker_id = $data[0];
		push @markerids, $marker_id;
		
		my @alleles = split('/', $data[1]);
		my $major_allele = $alleles[0]; #Major allele is ALWAYS leftmost allele in HapMap file
		
		unless ($major_allele =~ m/[ACTG+-]/) { print "$major_allele\t$."; }
		my $position = $data[3];
		print $GM "$marker_id\t$i\t$position\n";
		
		@data = @data[11..$#data];
		
		for my $h (0..$#data) {
		
			if ($data[$h] =~ m/[N$major_allele]/ ) { push @{ $numerical{$taxa[$h]} }, '2'; } #Conservatively impute missing data to major allele; use 2 for major allele
			elsif ($data[$h] =~ m/[KMRSWY0]/ ) { push @{ $numerical{$taxa[$h]} }, '1'; } #Heterozygotes are 1
			elsif ($data[$h] =~ m/[ACTG+-]/ ) { push @{ $numerical{$taxa[$h]} }, '0'; } #Homozygotes that don't match major allele (i.e. minor alleles) are 0
		
		}
		
	}
	
	close $GM;
	
	my $markerlist = join("\t", @markerids);
	print $GD "taxa\t$markerlist\n";
	
	for my $g (0..$#taxa) {
	
		my $output = join("\t", $taxa[$g],@{$numerical{$taxa[$g]}});
		print $GD "$output\n";
	
	}

	close $GD;

}