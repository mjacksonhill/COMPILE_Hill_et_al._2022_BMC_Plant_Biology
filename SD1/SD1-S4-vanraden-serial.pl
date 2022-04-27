#generate kinship matrices from HapMap files
#files are generated in K_chr style; i.e. each chromosome's kinship matrix is computed from the marker data of the other 9

#R installation path specified on line 54
#chromosome number specified on lines 25, 31

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
use File::Basename;
use File::Copy;
use File::Path;
use Cwd qw(abs_path);
my $dirname = dirname(abs_path($0));
$dirname =~ s/\\/\//g;

for my $i (1..10) {

	mkdir "$i", 0755;
	
	my $counter = 1;

	for my $j (1..10) {
	
		unless ($i == $j) {
		
			copy("$j.GD.txt", "$i/tmp$counter.GD");
			copy("$j.GM.txt", "$i/tmp$counter.GM");
			$counter++;
			
		}
		
	}

	open my $rfile, '>', "$i/R.txt";
	print $rfile qq`library(compiler)
library(gplots)
setwd("$dirname")
source("gapit_functions.txt")
setwd("$dirname/$i")
myGAPIT <- GAPIT(file.GD="tmp",file.GM="tmp",file.Ext.GD="GD",file.Ext.GM="GM",Geno.View.output=FALSE,kinship.algorithm="VanRaden",file.from=1,file.to=9)
q()
`;	

	close $rfile;
	my $command = 'C:/"Program Files"/R/R-3.4.2/bin/x64/Rscript.exe'." $dirname/$i/R.txt > nul 2>\&1";
	system("$command");
	
	open my $infile, '<', "$i/GAPIT.Kin.VanRaden.csv";
	open my $outfile, '>', "$i"."k.txt";
	
	while (<$infile>) {
	
		$_ =~ s/,/\t/g;
		print $outfile $_;
	
	}
	
	close $infile;
	close $outfile;
	rmtree("$i");
	
	print "Chromosome $i complete!\n";

}

while (1) {
  my $child = waitpid(-1, 0);
  last if $child == -1;
}
