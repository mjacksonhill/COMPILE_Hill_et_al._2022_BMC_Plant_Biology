#convert all VCF files in script folder to HapMap files

#lines 78 and 79 refer to specific taxa for these studies; to edit the entire VCF file without changing the contents:
	#change the range of positions in line 77 to read (0..n-1), where n is the total number of taxa
	#change the list of names in line 78 to a double-quoted, comma-separated, ordered list of the taxa in the VCF file
#total number of files to convert is specified on line 63 (e.g. for one VCF file per chromosome)
#line 118 filters marker data based on missing and minor allele proportions; edit or comment out this line as needed

#STRUCTURE
#SD1/
#├── SD1-S2.pl #script to perform renaming
#├── Input .vcf file(s) named 1.vcf, 2.vcf, etc.
#└── Output .hmp.txt file(s)

use strict;
use warnings;

my %IUPAC = (
	'AC' => 'M',
	'CA' => 'M',
	'AG' => 'R',
	'GA' => 'R',	
	'AT' => 'W',
	'TA' => 'W',
	'CG' => 'S',
	'GC' => 'S',	
	'CT' => 'Y',
	'TC' => 'Y',	
	'GT' => 'K',
	'TG' => 'K',	
	'-A' => '0',
	'A-' => '0',	
	'-C' => '0',
	'C-' => '0',	
	'-G' => '0',
	'G-' => '0',	
	'-T' => '0',
	'T-' => '0',	
	'+A' => '0',
	'A+' => '0',	
	'+C' => '0',
	'C+' => '0',	
	'+G' => '0',
	'G+' => '0',	
	'+T' => '0',
	'T+' => '0',	
	'+-' => '0',	
	'-+' => '0',	
	'AN' => 'N',
	'NA' => 'N',	
	'CN' => 'N',
	'NC' => 'N',	
	'GN' => 'N',
	'NG' => 'N',	
	'TN' => 'N',
	'NT' => 'N',	
	'+N' => 'N',
	'N+' => 'N',
	'-N' => 'N',
	'N-' => 'N',	
);

for my $i (1..10) {

open my $in, '<', "$i.vcf";
open my $out, '>', "$i.hmp.txt";

while (<$in>) {
	
	if ($_ =~ m/^#[^#]/) {
	
		last;
		
	}
	
}

my @positions = (23,916,917,918,919,920,921,922,123,923,924,925,926,927,928,929,930,931,932,933,934,935,936,937,938,939,940,941,942,943,944,945,946,947,948,949,950,951,952,953,954,955,956,957,958,146,959,960,961,962,963,971,964,965,966,967,968,969,970,972,973,974,975,976,977,212,978,979,980,981,982,983,218,984,985,986,987,988,989,990,991,992,993,994,995,996,997,998,999,1000,1001,224,1002,1003,1004,1005,1006,228,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,1034,1035,1036,1037,1049,1038,1040,1041,1042,1043,1044,1045,1046,1047,1048,1039,1050,1051,1052,1053,1054,1055,1056,1057,407,1058,1059,1060,1061,1062,1063,1064,430,461,462,1065,1066,1067,1068,1069,1070,1071,1072,1073,1074,1075,1076,1077,1078,1079,1080,517,1081,1082,1083,523,1084,1085,1086,1087,1088,1089,1090,1091,1092,1093,1094,1095,1096,1097,1098,1099,1100,1101,1102,1103,1104,1105,1106,1107,1108,1109,1110,1111,1112,1113,1114,1115,1116,1117,1118,1119,1120,1121,1122,1123,1124,1125,1126,1127,1128,1129,541,1130,1131,1132,1133,1134,1135,1136,1137,1138,1139,1140,1141,1142,1143,1144,1145,1146,1147,1148,1149,1150,1151,1152,1153,1154,1155,1156,1157,1158,1159,1160,1161,1162,1163,1164,1165,1166,1167,1168,1169,1170,1171,1172,1173,1174,1175,1176,725,732,1177,1178);
my @reference_taxa = ("207","33-16","38-11","4226","4722","A188","A214N","A239","A272","A441-5","A554","A556","A6","A619","A632","A634","A635","A641","A654","A659","A661","A679","A680","A682","Ab28A","B10","B103","B104","B105","B109","B14A","B164","B2","B37","B46","B52","B57","B64","B68","B73","B75","B76","B77","B79","B84","B97","C103","C123","C49A","CH701-30","CH9","CI.7","CI187-2","CI21E","CI28A","CI31A","CI3A","CI64","CI66","CI90C","CI91B","CM105","CM37","CM7","CML10","CML103","CML108","CML11","CML14","CML154Q","CML157Q","CML158Q","CML206","CML218","CML220","CML228","CML238","CML247","CML254","CML258","CML261","CML264","CML277","CML281","CML287","CML311","CML314","CML321","CML322","CML323","CML328","CML330","CML331","CML332","CML333","CML341","CML38","CML418","CML45","CML5","CML52","CML61","CML69","CML77","CML91","CML92","CMV3","CO106","CO125","CO255","D940Y","DE_2","DE1","DE811","E2558W","EP1","F2834T","F6","F7","GA209","GT112","H105W","H49","H84","H91","H95","H99","Hi27","HP301","Hy","I137TN","I205","I29","IA2132","Ia5125","IDS28","IDS69","IDS91","Il101","Il14H","Il677a","K148","K4","K55","K64","Ki11","Ki14","Ki2021","Ki21","Ki3","Ki43","Ki44","Ky21","KY226","KY228","L317","L578","LH132","LH74","LH82","M14","M162W","M37W","MEF156-55-2","Mo17","Mo18W","MO1W","Mo24W","Mo44","Mo45","Mo46","Mo47","MoG","Mp339","MS1334","MS153","MS71","Mt42","N192","N28Ht","N6","NC222","NC230","NC232","NC236","NC238","NC250","NC258","NC260","NC262","NC264","NC290A","NC294","NC296","NC296A","NC298","NC300","NC302","NC304","NC306","NC310","NC314","NC318","NC320","NC324","NC326","NC328","NC33","NC336","NC338","NC340","NC342","NC344","NC346","NC348","NC350","NC352","NC354","NC356","NC358","NC360","NC362","NC364","NC366","NC368","ND246","Oh40B","Oh43","Oh43E","Oh603","OH7B","Os420","P39","Pa762","Pa875","Pa880","Pa91","R168","R177","R229","R4","SA24","SC213R","SC357","SC55","SD40","SD44","Sg1533","Sg18","T232","T234","T8","Tx303","Tx601","Tzi10","Tzi11","Tzi16","Tzi18","Tzi25","Tzi8","Tzi9","U267Y","VA102","Va14","Va17","Va22","Va26","Va35","Va59","Va85","Va99","VaW6","W117Ht","W153R","W182B","W22","W64A","Wf9","Yu796_NS");

my $hmp_taxa = join("\t", @reference_taxa);
print $out "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\t$hmp_taxa\n";

while (<$in>) {

	chomp $_;
	my @fields = split("\t", $_);
	
	my $chrom = $fields[0];
	my $pos = $fields[1];
	my $id = $fields[2];
	my $ref_allele = $fields[3];
	
	my $alt_allele_field = $fields[4];
	$alt_allele_field =~ s/<INS>/\+/;
	$alt_allele_field =~ s/<DEL>/-/;
	my @alt_alleles = split(',', $alt_allele_field);
	my $alt1 = $alt_alleles[0];
	my $alt2 = $alt_alleles[1] // undef;
	my $alt3 = $alt_alleles[2] // undef;
	my $alleles = join('/', $ref_allele, @alt_alleles);
	
	@fields = @fields[9..$#fields];
	
	my @data;
	for (@positions) {
	
		push (@data, $fields[$_]);
	
	}
	
	my $data_field = join("\t", @data[0..$#data]);
	$data_field =~ s|\./\.|N|g;
	
	my $missing = () = $data_field =~ m|\tN\t|g;
	my $minor = () = $data_field =~ m`(/1|1/|/2|2/|3/|/3)`g;
	
	if (($missing > 69) || (($minor/(558-(2*$missing)) < 0.05))) { next; } 
	
	$data_field =~ s`:.*?(\t|\n)`$1`g;
	$data_field =~ s|:.*?$||g;
	$data_field =~ s|0/0|$ref_allele|g;
	$data_field =~ s|1/1|$alt1|g;
	
	my $zeroone = $IUPAC{"$ref_allele"."$alt1"};
	$data_field =~ s`(0/1|1/0)`$zeroone`g;
	
	if (defined($alt2)) { 
	
		$data_field =~ s|2/2|$alt2|g; 
		my $zerotwo = $IUPAC{"$ref_allele"."$alt2"};
		my $onetwo = $IUPAC{"$alt1"."$alt2"};
		$data_field =~ s`(0/2|2/0)`$zerotwo`g;
		$data_field =~ s`(1/2|2/1)`$onetwo`g;
		
	}
	
	if (defined($alt3)) { 
	
		$data_field =~ s|3/3|$alt3|g; 
		my $zerothree = $IUPAC{"$ref_allele"."$alt3"};
		my $onethree = $IUPAC{"$alt1"."$alt3"};
		my $threetwo = $IUPAC{"$alt3"."$alt2"};
		$data_field =~ s`(0/3|3/0)`$zerothree`g;
		$data_field =~ s`(1/3|3/1)`$onethree`g;
		$data_field =~ s`(3/2|2/3)`$threetwo`g;
	
	}
	
	print $out "$id\t$alleles\t$chrom\t$pos\t+\tNA\tNA\tNA\tNA\tNA\tNA\t$data_field\n";
	
}

close $in;
close $out;

}