use strict;
use warnings;
use File::Basename;


my ($bam,$JX) = @ARGV;

my $cutoff = 1000; # 1000bp

my ($del_sp,$del_ep); # start/end
if ($JX eq "JX_2"){
  $del_sp = 34232240;
  $del_ep = 34241309;
  #$region = "chr18:34232240-34241309";
}else{
  $del_sp = 34255473;
  $del_ep = 34262838;
  #$region = "chr18:34255473-34262838";
}


my $samtools = "samtools";
open IN, "$samtools view $bam |" or die;
while (<IN>){
  chomp;
  my @arr = split /\t/, $_;
  my $seq_name = $arr[0];
  my $flag = $arr[1];
  my $pos = $arr[3];
  my $MAPQ = $arr[4]; # default 60
  my $cigar = $arr[5];
  my $mate_pos = $arr[7];
  my $insert_size = $arr[8];

  # 排除 secondary / supplementary
  next if ($flag & 0x100);
  next if ($flag & 0x800);

  # 判断：起点或终点落在del范围1000bp以内
  if ($pos < $mate_pos){
    if (abs($pos-$del_sp) <= $cutoff and abs($mate_pos-$del_ep) <= $cutoff){
      # OK
      print "$_\n";
    }
  }else{
    if (abs($pos-$del_ep) <= $cutoff and abs($mate_pos-$del_sp) <= $cutoff){
      # OK
      print "$_\n";
    }
  }

}

close IN;


