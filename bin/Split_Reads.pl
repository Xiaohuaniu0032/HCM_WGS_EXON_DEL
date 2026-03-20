use strict;
use warnings;
use Getopt::Long;

# ------------------------------------------------------------
#   - Read BAM via samtools view
#   - Keep primary alignments only (not secondary 0x100, not supplementary 0x800)
#   - Require READ1 or READ2 and SA tag
#   - Parse SA CIGAR
# ------------------------------------------------------------

my ($bam,$JX) = @ARGV;

my $cutoff = 1; # 1bp => breakpoint is base level

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


# store PE read info
my %info;
my @soft_clip_reads;

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

  push @{$info{$seq_name}}, $_;

  # soft cliped reads / NOT SPLIT READS
  if ($cigar =~ /(\d+)M(\d+)S/){
    # left breakpoint
    my $end_pos = $pos + $1 - 1;
    if ($end_pos + 1 == $del_sp){
      # OK
      push @soft_clip_reads, $seq_name;
    }
  }elsif ($cigar =~ /(\d+)S(\d+)M/){
    # right breakpoint
    my $end_pos = $pos;
    if ($end_pos - 1 == $del_ep){
      # OK
      push @soft_clip_reads, $seq_name;
    }
  }else{
    next;
  }
}
close IN;


for my $seq_name (@soft_clip_reads){
  my $READ1 = $info{$seq_name}->[0];
  my $READ2 = $info{$seq_name}->[1];
  if ($READ1){
    print "$READ1\n";
  }
  if ($READ2){
    print "$READ2\n";
  }
}
