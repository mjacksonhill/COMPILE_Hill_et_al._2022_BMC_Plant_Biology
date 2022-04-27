#standardize maize phenotype data names to those of Romay et al., 2013

#STRUCTURE
#SD1/
#├── SD1-S1.pl #script to perform renaming
#├── Input X.txt file(s) #phenotype data with tab-separated data with genotype names in first column
#├── Output X-named.txt file(s) #renamed phenotype data
#└── SD1-Files/
#    └── Romay.txt #contains standardized line names according to Romay et al. 2013

use strict;
use warnings;

open my $Romay, '<', 'SD1-Files/Romay.txt';
chomp(my @Romay = <$Romay>);
close $Romay;
my $total = $#Romay+1;

my @files = glob "*.txt";
for my $file (@files) {

	my %queries;

	open my $in, '<', "$file";
	my $header = <$in>;
	
	while (<$in>) {
	
		#certain renames specific to these studies
		chomp $_;
		my @terms = split("\t", $_);
		$terms[0] =~ s/:.*?$//;
		$terms[0] =~ s/B2-good/B2/;
		$terms[0] =~ s/Il677A/Il677a/;
		$terms[0] =~ s/Ky226/KY226/;
		$terms[0] =~ s/Ky228/KY228/;
		$terms[0] =~ s/MO17/Mo17/;
		$terms[0] =~ s/Mo1W/MO1W/;
		$terms[0] =~ s/OH43/Oh43/;
		$terms[0] =~ s/PA91/Pa91/;
		$terms[0] =~ s/TX601/Tx601/;
		$terms[0] =~ s/TzI10/Tzi10/;
		$terms[0] =~ s/TzI18/Tzi18/;
		$terms[0] =~ s/TzI8/Tzi8/;
		$terms[0] =~ s/Va102/VA102/;
		$terms[0] =~ s/VA26/Va26/;
		$terms[0] =~ s/W117HT/W117Ht/;
		$terms[0] =~ s/WF9/Wf9/;
		$terms[0] =~ s/W22R-r-rstd/W22R-r-std/;
		$terms[0] =~ s/W22R-r-std_CS-2909-1/W22R-r-std/;
		$terms[0] =~ s/trisacum/tripsacum/;
		$queries{$terms[0]} = $terms[1];
	
	}
	
	close $in;

	open my $out, '>', "$file-named.txt";
	print $out $header;
	
	my $i = 0;
	
	for my $reference (@Romay) {
	
		my $data = $queries{$reference} // "NaN";
		if ($data eq "NaN") { $i++; }
		print $out "$reference\t$data\n";
	
	}

	close $out;
	print "Missing taxa were $i of $total for $file\n"

}