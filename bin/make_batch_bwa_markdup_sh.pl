#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);
use File::Path qw(make_path);

# ============================================================
# Script:
#   make_batch_bwa_markdup_sh.pl
#
# Function:
#   Scan subdirectories under a FASTQ simulation directory.
#   For each sample directory, detect paired FASTQ files automatically,
#   create an output subdirectory named by sample,
#   and generate one BWA + samtools markdup shell script.
#
# Expected input directory structure:
#   random_gene_del_sim/
#     HCM_FHOD3_DEL_5000_10000_rep2/
#       HCM_FHOD3_DEL_5000_10000_rep2.R1.fastq
#       HCM_FHOD3_DEL_5000_10000_rep2.R2.fastq
#       HCM_FHOD3_DEL_5000_10000_rep2.truth.txt
#
# Output structure:
#   outdir/
#     HCM_FHOD3_DEL_5000_10000_rep2/
#       HCM_FHOD3_DEL_5000_10000_rep2.bwa_markdup.sh
#       run after sh:
#         HCM_FHOD3_DEL_5000_10000_rep2.sorted.markdup.bam
#         HCM_FHOD3_DEL_5000_10000_rep2.sorted.markdup.bam.bai
#         HCM_FHOD3_DEL_5000_10000_rep2.markdup.metrics.txt
#         HCM_FHOD3_DEL_5000_10000_rep2.bwa_markdup.log
#
# Example:
#   perl make_batch_bwa_markdup_sh.pl \
#     --indir /ehpcdata/fulongfei/project/SV_Caller_Test_20260522/simulate_fq/random_gene_del_sim \
#     --ref /path/to/hg19.fa \
#     --outdir /ehpcdata/fulongfei/project/SV_Caller_Test_20260522/data/BAM/random_gene_del_sim \
#     --threads 8 \
#     --bwa /path/to/bwa \
#     --samtools /path/to/samtools
# ============================================================

my ($indir, $ref, $outdir);
my $threads  = 8;
my $bwa      = "bwa";
my $samtools = "samtools";
my $platform = "ILLUMINA";
my $keep_tmp = 0;
my $force    = 0;
my $run_list = "";
my $help     = 0;

GetOptions(
    "indir=s"    => \$indir,
    "ref=s"      => \$ref,
    "outdir=s"   => \$outdir,
    "threads=i"  => \$threads,
    "bwa=s"      => \$bwa,
    "samtools=s" => \$samtools,
    "platform=s" => \$platform,
    "keep-tmp!"  => \$keep_tmp,
    "force!"     => \$force,
    "run-list=s" => \$run_list,
    "help"       => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless $indir && $ref && $outdir;

check_dir($indir, "Input FASTQ root directory");
check_file($ref, "Reference FASTA");

$indir = abs_path($indir);
$ref   = abs_path($ref);

make_path($outdir) unless -d $outdir;
$outdir = abs_path($outdir);

$run_list ||= "$outdir/run_all_bwa_markdup.sh";

open my $RUN, ">", $run_list or die "[ERROR] Cannot write run list: $run_list\n";
print $RUN "#!/usr/bin/env bash\n";
print $RUN "set -euo pipefail\n\n";

my @sample_dirs = get_sample_dirs($indir);

die "[ERROR] No sample subdirectories found under: $indir\n" unless @sample_dirs;

my $generated = 0;
my $skipped   = 0;

foreach my $sample_dir (@sample_dirs) {

    my $sample = basename($sample_dir);

    my ($fq1, $fq2) = find_pair_fastq($sample_dir);

    unless ($fq1 && $fq2) {
        warn "[WARN] Skip sample because paired FASTQ not found: $sample_dir\n";
        $skipped++;
        next;
    }

    my $sample_outdir = "$outdir/$sample";
    make_path($sample_outdir) unless -d $sample_outdir;

    my $sh = "$sample_outdir/$sample.bwa_markdup.sh";

    if (-e $sh && !$force) {
        warn "[WARN] Shell script exists, skip: $sh\n";
        $skipped++;
        next;
    }

    write_bwa_markdup_sh(
        fq1        => $fq1,
        fq2        => $fq2,
        ref        => $ref,
        sample     => $sample,
        outdir     => $sample_outdir,
        sh         => $sh,
        threads    => $threads,
        bwa        => $bwa,
        samtools   => $samtools,
        platform   => $platform,
        keep_tmp   => $keep_tmp,
    );

    print $RUN "bash '$sh'\n";
    $generated++;

    print "[INFO] Generated: $sh\n";
    print "[INFO]   Sample: $sample\n";
    print "[INFO]   FQ1   : $fq1\n";
    print "[INFO]   FQ2   : $fq2\n";
}

close $RUN;

chmod 0755, $run_list or die "[ERROR] Cannot chmod run list: $run_list\n";

print "\n[INFO] Batch shell generation finished.\n";
print "[INFO] Input root      : $indir\n";
print "[INFO] Output root     : $outdir\n";
print "[INFO] Generated       : $generated\n";
print "[INFO] Skipped         : $skipped\n";
print "[INFO] Run list        : $run_list\n";
print "[INFO] Run all command : bash $run_list\n";


# ------------------------------------------------------------
# Subroutines
# ------------------------------------------------------------

sub get_sample_dirs {
    my ($root) = @_;

    opendir my $DH, $root or die "[ERROR] Cannot open directory: $root\n";

    my @dirs;
    while (my $item = readdir $DH) {
        next if $item eq "." || $item eq "..";

        my $path = "$root/$item";
        next unless -d $path;

        push @dirs, abs_path($path);
    }

    closedir $DH;

    @dirs = sort @dirs;
    return @dirs;
}


sub find_pair_fastq {
    my ($dir) = @_;

    opendir my $DH, $dir or die "[ERROR] Cannot open sample directory: $dir\n";

    my @files;
    while (my $item = readdir $DH) {
        next if $item eq "." || $item eq "..";

        my $path = "$dir/$item";
        next unless -f $path;

        next unless $item =~ /\.(fastq|fq)(\.gz)?$/i;

        push @files, abs_path($path);
    }

    closedir $DH;

    @files = sort @files;

    my @r1;
    my @r2;

    foreach my $file (@files) {
        my $base = basename($file);

        if ($base =~ /(^|[._-])R?1([._-]|$)/i) {
            push @r1, $file;
        }
        elsif ($base =~ /(^|[._-])R?2([._-]|$)/i) {
            push @r2, $file;
        }
    }

    if (@r1 == 1 && @r2 == 1) {
        return ($r1[0], $r2[0]);
    }

    if (@r1 > 1 || @r2 > 1) {
        warn "[WARN] Multiple R1/R2 FASTQ files detected in $dir\n";
        warn "[WARN] R1 files:\n";
        warn "       $_\n" for @r1;
        warn "[WARN] R2 files:\n";
        warn "       $_\n" for @r2;
        return;
    }

    return;
}


sub write_bwa_markdup_sh {
    my (%args) = @_;

    my $fq1      = $args{fq1};
    my $fq2      = $args{fq2};
    my $ref      = $args{ref};
    my $sample   = $args{sample};
    my $outdir   = $args{outdir};
    my $sh       = $args{sh};
    my $threads  = $args{threads};
    my $bwa      = $args{bwa};
    my $samtools = $args{samtools};
    my $platform = $args{platform};
    my $keep_tmp = $args{keep_tmp};

    check_file($fq1, "FASTQ R1");
    check_file($fq2, "FASTQ R2");

    my $tmp_dir          = "$outdir/tmp";
    my $name_sorted_bam  = "$tmp_dir/$sample.name_sorted.bam";
    my $fixmate_bam      = "$tmp_dir/$sample.fixmate.bam";
    my $coord_sorted_bam = "$tmp_dir/$sample.coord_sorted.bam";
    my $markdup_bam      = "$outdir/$sample.sorted.markdup.bam";
    my $markdup_metrics  = "$outdir/$sample.markdup.metrics.txt";
    my $log_file         = "$outdir/$sample.bwa_markdup.log";

    my $rg = build_read_group(
        sample   => $sample,
        platform => $platform,
        library  => $sample,
        unit     => $sample,
    );

    open my $SH, ">", $sh or die "[ERROR] Cannot write shell script: $sh\n";

    print $SH <<"SH";
#!/usr/bin/env bash
set -euo pipefail

# Auto-generated by make_batch_bwa_markdup_sh.pl

FQ1='$fq1'
FQ2='$fq2'
REF='$ref'
SAMPLE='$sample'
OUTDIR='$outdir'
TMPDIR='$tmp_dir'
THREADS='$threads'

BWA='$bwa'
SAMTOOLS='$samtools'

NAME_SORTED_BAM='$name_sorted_bam'
FIXMATE_BAM='$fixmate_bam'
COORD_SORTED_BAM='$coord_sorted_bam'
MARKDUP_BAM='$markdup_bam'
MARKDUP_METRICS='$markdup_metrics'
LOG='$log_file'

mkdir -p "\$OUTDIR"
mkdir -p "\$TMPDIR"

echo "[INFO] Sample: \$SAMPLE" > "\$LOG"
echo "[INFO] FQ1: \$FQ1" >> "\$LOG"
echo "[INFO] FQ2: \$FQ2" >> "\$LOG"
echo "[INFO] Reference: \$REF" >> "\$LOG"
echo "[INFO] Output BAM: \$MARKDUP_BAM" >> "\$LOG"
echo "[INFO] Threads: \$THREADS" >> "\$LOG"
echo "[INFO] Start time: \$(date)" >> "\$LOG"

echo "[INFO] Step 1: BWA alignment and name sorting" | tee -a "\$LOG"
"\$BWA" mem -t "\$THREADS" -R '$rg' "\$REF" "\$FQ1" "\$FQ2" 2>> "\$LOG" | \\
  "\$SAMTOOLS" sort -n -@ "\$THREADS" -o "\$NAME_SORTED_BAM" - 2>> "\$LOG"

echo "[INFO] Step 2: Fix mate information" | tee -a "\$LOG"
"\$SAMTOOLS" fixmate -@ "\$THREADS" -m "\$NAME_SORTED_BAM" "\$FIXMATE_BAM" 2>> "\$LOG"

echo "[INFO] Step 3: Coordinate sorting" | tee -a "\$LOG"
"\$SAMTOOLS" sort -@ "\$THREADS" -o "\$COORD_SORTED_BAM" "\$FIXMATE_BAM" 2>> "\$LOG"

echo "[INFO] Step 4: Mark duplicates" | tee -a "\$LOG"
"\$SAMTOOLS" markdup -@ "\$THREADS" -s "\$COORD_SORTED_BAM" "\$MARKDUP_BAM" > "\$MARKDUP_METRICS" 2>> "\$LOG"

echo "[INFO] Step 5: BAM indexing" | tee -a "\$LOG"
"\$SAMTOOLS" index -@ "\$THREADS" "\$MARKDUP_BAM" 2>> "\$LOG"

SH

    unless ($keep_tmp) {
        print $SH <<"SH";
echo "[INFO] Step 6: Remove temporary files" | tee -a "\$LOG"
rm -f "\$NAME_SORTED_BAM" "\$FIXMATE_BAM" "\$COORD_SORTED_BAM"
rmdir "\$TMPDIR" 2>/dev/null || true

SH
    }

    print $SH <<"SH";
echo "[INFO] Finished time: \$(date)" >> "\$LOG"
echo "[INFO] Done."
echo "[INFO] Final BAM: \$MARKDUP_BAM"
echo "[INFO] BAM index: \$MARKDUP_BAM.bai"
echo "[INFO] Markdup metrics: \$MARKDUP_METRICS"
echo "[INFO] Log file: \$LOG"

SH

    close $SH;

    chmod 0755, $sh or die "[ERROR] Cannot chmod shell script: $sh\n";
}


sub build_read_group {
    my (%args) = @_;

    my $sample   = $args{sample};
    my $platform = $args{platform} || "ILLUMINA";
    my $library  = $args{library}  || $sample;
    my $unit     = $args{unit}     || $sample;

    my $id = $unit;

    my $rg = "\@RG";
    $rg .= "\\tID:$id";
    $rg .= "\\tSM:$sample";
    $rg .= "\\tLB:$library";
    $rg .= "\\tPL:$platform";
    $rg .= "\\tPU:$unit";

    return $rg;
}


sub check_file {
    my ($file, $name) = @_;

    die "[ERROR] $name not found: $file\n" unless -e $file;
    die "[ERROR] $name is empty: $file\n" unless -s $file;
}


sub check_dir {
    my ($dir, $name) = @_;

    die "[ERROR] $name not found: $dir\n" unless -d $dir;
}


sub usage {
    return <<"USAGE";

Usage:
  perl make_batch_bwa_markdup_sh.pl \\
    --indir random_gene_del_sim \\
    --ref hg19.fa \\
    --outdir ./BAM

Required parameters:
  --indir      Root directory containing sample FASTQ subdirectories
  --ref        Reference genome FASTA
  --outdir     Output root directory

Optional parameters:
  --threads    Number of threads, default: 8
  --bwa        Path to bwa, default: bwa
  --samtools   Path to samtools, default: samtools
  --platform   Sequencing platform, default: ILLUMINA
  --keep-tmp   Keep temporary BAM files, default: off
  --force      Overwrite existing shell scripts, default: off
  --run-list   Output run-all shell script, default: <outdir>/run_all_bwa_markdup.sh
  --help       Show help information

Input example:
  random_gene_del_sim/
    HCM_FHOD3_DEL_5000_10000_rep2/
      HCM_FHOD3_DEL_5000_10000_rep2.R1.fastq
      HCM_FHOD3_DEL_5000_10000_rep2.R2.fastq
      HCM_FHOD3_DEL_5000_10000_rep2.truth.txt

Output example:
  outdir/
    HCM_FHOD3_DEL_5000_10000_rep2/
      HCM_FHOD3_DEL_5000_10000_rep2.bwa_markdup.sh

Run:
  bash <outdir>/run_all_bwa_markdup.sh

USAGE
}

