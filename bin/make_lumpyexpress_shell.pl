use strict;
use warnings;
use Cwd qw(abs_path);
use File::Basename qw(dirname basename);


my $binDIR = "/home/fulongfei/miniconda3/envs/lumpy_py27/bin";
my $lumpyexpress = "$binDIR/lumpyexpress";
my $config = "$binDIR/lumpyexpress.config";
my $samtools = "/opt/PUB/software/Bin/samtools-1.16.1/bin/samtools";

# param
my ($bam,$name,$outdir) = @ARGV;

my $runsh = "$outdir/run_lumpyexpress\_$name\.sh";
open SH, ">$runsh" or die;

# prepare discordant bam
my $discordants_bam = "$outdir/$name\.discordants.bam";
my $cmd = "$samtools view -b -F 1294 $bam \>$discordants_bam";
print SH "$cmd\n";

# prepare split bam
my $split_bam = "$outdir/$name\.split.bam";
$cmd = "$samtools view -h $bam | $binDIR/extractSplitReads_BwaMem -i stdin | samtools view -Sb - \>$split_bam";
print SH "$cmd\n";

# full bam has been sort by pos
# discordant/split bam do not need to sort again

# call lumpyexpress
my $vcf = "$outdir/$name\.lumpy.vcf";
$cmd = "$lumpyexpress -B $bam -S $split_bam -D $discordants_bam -o $vcf";
print SH "$cmd\n";
close SH;
`chmod 755 $runsh`;

