#generate kinship matrices from HapMap files
#files are generated in K_chr style; i.e. each chromosome's kinship matrix is computed from the marker data of the other 9

#chromosome number specified on line 18

#STRUCTURE
#SD1/ #vanraden scripts use vanraden method and require R and relevant packages
#├── SD1-S4.pl #script to perform parallel calculation by loiselle method
#├── SD1-S4-serial.pl #script to perform serial calculation by loiselle method (less computationally expensive, much longer runtime)
#├── SD1-S4-vanraden.pl #script to perform parallel calculation by vanraden (overall much less expensive)
#├── SD1-S4-vanraden-serial.pl #script to perform serial calculation by vanraden (yet less expensive, still longer runtime)
#├── Input 1.hmp.txt - ##.hmp.txt genotype files
#└── Output 1.k.txt - ##.k.txt kinship matrix files

use strict;
use warnings;

my $chromosome_number = 10;
my $debug = 0;

my %IUPAC = ( #For converting single nucleotide IUPAC codes
	'A' => ['A','A'],
	'C' => ['C','C'],
	'G' => ['G','G'],
	'T' => ['T','T'],
	'+' => ['+','+'],
	'-' => ['-','-'],
	'K' => ['G','T'],
	'M' => ['A','C'],
	'R' => ['A','G'],
	'S' => ['C','G'],
	'W' => ['A','T'],
	'Y' => ['C','T'],
	'0' => ['+','-'],
	'N' => ['N','N'],
);

sub iupac {
	
   my $letter_1 = @{ $IUPAC{$_[0]} }[0];
   my $letter_2 = @{ $IUPAC{$_[0]} }[1];
   return ("$letter_1","$letter_2")
   
}

open my $single_denoms, '>', "single_denominators.txt";
my @taxa_names;

for my $i (1..$chromosome_number) {

die "could not fork" unless defined(my $prepid = fork); #Parallel
unless ($prepid) {

	my %store;
	my %kinship_numerators;
	open my $infile, '<', "$i.hmp.txt";
	my $header = <$infile>;
	
	my $denominator = 0;
	
	while (<$infile>) {
		my @terms = split("\t", $_);
		chomp @terms;
		
		my @marker_data = @terms[11..$#terms];
		my @marker_data_pruned = grep(!/N/,@marker_data);
		my $total_alleles = ($#marker_data_pruned+1)*2;
		
		my %allele_numbers; #For storing numbers of alleles

		for my $marker (@marker_data_pruned) { #Tally up single-letter allele counts from HapMap row
		
			my @letters = iupac($marker);
			$allele_numbers{$letters[0]}++; 
			$allele_numbers{$letters[1]}++; 
			
		}
		
		my $fA = ($allele_numbers{'A'} // 0)/$total_alleles; #Final allele frequency terms
		my $fC = ($allele_numbers{'C'} // 0)/$total_alleles;
		my $fG = ($allele_numbers{'G'} // 0)/$total_alleles;
		my $fT = ($allele_numbers{'T'} // 0)/$total_alleles;
		my $fPlus = ($allele_numbers{'+'} // 0)/$total_alleles;
		my $fMinus = ($allele_numbers{'-'} // 0)/$total_alleles;
		
		my $term_1a = ($fA*(1-$fA))+($fC*(1-$fC))+($fG*(1-$fG))+($fT*(1-$fT))+($fPlus*(1-$fPlus))+($fMinus*(1-$fMinus));
		my $term_1 = $term_1a/($total_alleles-1);
		$denominator += $term_1a;
		
		push @{ $store{'A'} }, "$fA";
		push @{ $store{'C'} }, "$fC";
		push @{ $store{'G'} }, "$fG";
		push @{ $store{'T'} }, "$fT";
		push @{ $store{'Plus'} }, "$fPlus";
		push @{ $store{'Minus'} }, "$fMinus";
		push @{ $store{'Term1'} }, "$term_1";
		
	}
	
	close $infile;
	print $single_denoms "$i\t$denominator\n";
	
	my %hapmap;
	open my $hmp, '<', "$i.hmp.txt";
	my $hmphead = <$hmp>;
	my @taxa = split("\t",$hmphead);
	@taxa = @taxa[11..$#taxa];
	chomp @taxa;
	@taxa_names = @taxa;
	
	while (<$hmp>) {
	
		my @SNPs = split("\t",$_);
		@SNPs = @SNPs[11..$#SNPs];
		chomp @SNPs;
	
		for my $m (0..$#SNPs) {
		
			push @{ $hapmap{$taxa[$m]} }, $SNPs[$m];

		}
	
	}
	
	close $hmp;
	
	my $total_ops = (($#taxa*($#taxa+1))/2);
	my $counter = 1;
	
	for my $j (0..$#taxa) { #For every taxon
	
		for my $k ($j..$#taxa) { #For every taxon we haven't already looped through
			
			my @j_markers = @{ $hapmap{$taxa[$j]} };
			chomp @j_markers;
			my @k_markers = @{ $hapmap{$taxa[$k]} };
			chomp @k_markers;
			
			my $Term2 = 0;
		
			for my $m (0..$#j_markers) { #Across all markers
				
				my $pA = @{ $store{'A'} }[$m]; #Grab allele frequencies from earlier step
				my $pC = @{ $store{'C'} }[$m];
				my $pG = @{ $store{'G'} }[$m];
				my $pT = @{ $store{'T'} }[$m];
				my $pPlus = @{ $store{'Plus'} }[$m];
				my $pMinus = @{ $store{'Minus'} }[$m];
				my $Term1 = @{ $store{'Term1'} }[$m];
				
				my %m; #Populates hash with allele numbers for both taxa
				my @j = iupac($j_markers[$m]);
				my @k = iupac($k_markers[$m]);
				$m{'j'}{$j[0]} += 0.5;
				$m{'k'}{$k[0]} += 0.5;
				$m{'j'}{$j[1]} += 0.5;
				$m{'k'}{$k[1]} += 0.5;
				
				my $jA = $m{'j'}{'A'} // 0;
				my $jC = $m{'j'}{'C'} // 0;
				my $jG = $m{'j'}{'G'} // 0;
				my $jT = $m{'j'}{'T'} // 0;
				my $jPlus = $m{'j'}{'+'} // 0;
				my $jMinus = $m{'j'}{'-'} // 0;
				
				my $kA = $m{'k'}{'A'} // 0;
				my $kC = $m{'k'}{'C'} // 0;
				my $kG = $m{'k'}{'G'} // 0;
				my $kT = $m{'k'}{'T'} // 0;
				my $kPlus = $m{'k'}{'+'} // 0;
				my $kMinus = $m{'k'}{'-'} // 0;
				
				$Term2 = $Term2+(((($jA-$pA)*($kA-$pA))+(($jC-$pC)*($kC-$pC))+(($jG-$pG)*($kG-$pG))+(($jT-$pT)*($kT-$pT))+(($jPlus-$pPlus)*($kPlus-$pPlus))+(($jMinus-$pMinus)*($kMinus-$pMinus)))+$Term1);
				$kinship_numerators{$taxa[$j]}{$taxa[$k]} = $Term2; #Square matrix in hash
				$kinship_numerators{$taxa[$k]}{$taxa[$j]} = $Term2;
				
			}
			
			$counter++;
		
		}
		
	}
	
	open my $numerator_file, '>', "$i.numerator.txt";
	my @list = sort keys(%kinship_numerators);
	
	for my $label (@list) {
	
		print $numerator_file "$label";
		
		for my $u (0..$#list) {
			
			print $numerator_file "\t$kinship_numerators{$label}{$list[$u]}";

		}
		
		print $numerator_file "\n";
		
	}
	
	close $numerator_file;
	exit;

}

}

while (1) {
  my $child = waitpid(-1, 0);
  last if $child == -1;
}

close $single_denoms;
my %kchr_denominators;

my $total_denom = 0;
my %single_denominators;
open my $singles, '<', "single_denominators.txt";

while (<$singles>) {

	chomp $_;
	my @terms = split("\t", $_);
	$single_denominators{$terms[0]} = $terms[1];

}

close $singles;

for my $p (1..$chromosome_number) {

    $total_denom += $single_denominators{$p};
	
}

print "Performing K_Chr Arithmetic...\n";

for my $j (1..10) { #Print k_chr denominators-- sum of denominators for all chromosomes except the one in question

	$kchr_denominators{$j} = $total_denom-$single_denominators{$j};
	
}

my %kchr; #initialize matrix for adding up numerators

for my $name (@taxa_names) {

	for my $num (0..$#taxa_names) {
	
		$kchr{$name}[$num] = 0;
	
	}

}

for my $h (1..$chromosome_number) {

	open my $infile, '<', "$h.numerator.txt";
	
	while (<$infile>) {
	
		my @data = split("\t", $_);
		chomp @data;
		my $name = shift(@data);
		
		for my $col (0..$#data) {
		
			$kchr{$name}[$col] += $data[$col];
		
		}
		
	}

}

print "Exporting Final Files...\n";

for my $v (1..$chromosome_number) {

	open my $infile, '<', "$v.numerator.txt";
	open my $outfile, '>', "$v"."k.txt";
	
	while (<$infile>) {
	
		my @data = split("\t", $_);
		chomp @data;
		my $name = shift(@data);
		print $outfile "$name";
		
		for my $q (0..$#data) {
		
			my $total = $kchr{$name}[$q];
			my $minus = $data[$q];
			my $value = ($total-$minus)/$kchr_denominators{$v};
			$value = sprintf("%.9f", $value); #Round to 9 decimal places
			
			print $outfile "\t$value";
		
		}
		
		print $outfile "\n";
	
	}
	
	close $outfile;

}

unless ($debug == 1) {

	unlink "single_denominators.txt";

	for my $t (1..$chromosome_number) {

		unlink "$t.numerator.txt";

	}

}