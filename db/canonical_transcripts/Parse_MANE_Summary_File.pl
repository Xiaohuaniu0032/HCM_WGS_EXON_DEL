use strict;
use warnings;
use File::Basename;


my $file = "MANE.GRCh38.v1.4.summary.txt";

print "name\tRefSeq_prot\tMANE_status\n";

open IN, "$file" or die;
<IN>;
while (<IN>){
	chomp;
	my @arr = split /\t/, $_;
	my $name = $arr[3];
	my $RefSeq_prot = $arr[5]; # NM_130786.4
	$RefSeq_prot =~ s/\.(\d+)$//;
	my $MANE_status = $arr[-5]; # MANE Select / MANE Plus Clinical
	print "$name\t$RefSeq_prot\t$MANE_status\n";
}
close IN;



