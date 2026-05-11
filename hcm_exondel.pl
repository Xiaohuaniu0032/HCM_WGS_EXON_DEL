#!/usr/bin/env perl
# -*- coding: utf-8 -*-

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(dirname basename);
use File::Path qw(make_path);
use Cwd qw(abs_path);
use FindBin qw($Bin);

# ============================================================
# HCMExonDel main pipeline
# Author: Fulongfei Fu
# Description:
#   Main control script for detecting exon-level deletions in
#   HCM-related genes from coordinate-sorted WGS BAM files.
#
# Core evidence:
#   1. Gene/exon depth reduction
#   2. Discordant read pairs
#   3. Split reads / soft-clipped reads
# ============================================================

my $config;
my $bam;
my $bam_list;
my $sample;
my $outdir = "results";
my $threads;
my $target_bed;
my $force = 0;
my $keep_tmp = 0;
my $help = 0;

GetOptions(
    "config=s"     => \$config,
    "bam=s"        => \$bam,
    "bam-list=s"   => \$bam_list,
    "sample=s"     => \$sample,
    "outdir=s"     => \$outdir,
    "threads=i"    => \$threads,
    "target-bed=s" => \$target_bed,
    "force"        => \$force,
    "keep-tmp"     => \$keep_tmp,
    "help"         => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

# -----------------------------
# Check required arguments
# -----------------------------
die usage() unless $config;
die "[ERROR] Please provide --bam or --bam-list\n" unless ($bam || $bam_list);
die "[ERROR] Please use either --bam or --bam-list, not both\n" if ($bam && $bam_list);
die "[ERROR] --sample is required when using --bam\n" if ($bam && !$sample);

$config = abs_path($config);
die "[ERROR] Config file not found: $config\n" unless -s $config;

# -----------------------------
# Read config file
# -----------------------------
my %CONF = read_config($config);

# Command-line parameters have higher priority than config file
$threads    ||= $CONF{THREADS} || 4;
$target_bed ||= $CONF{TARGET_REGION_BED} || $CONF{EXON_ANNOTATION_BED};
$keep_tmp   ||= $CONF{KEEP_TMP} || 0;

die "[ERROR] TARGET_REGION_BED or EXON_ANNOTATION_BED is required in config or --target-bed\n"
    unless $target_bed;

$target_bed = abs_path($target_bed);
die "[ERROR] Target BED file not found: $target_bed\n" unless -s $target_bed;

$outdir = abs_path($outdir) || $outdir;
make_path($outdir) unless -d $outdir;

# -----------------------------
# Check dependent scripts
# -----------------------------
my $script_depth      = "$Bin/run_gene_mean_depth.pl";
my $script_discordant = "$Bin/run_discordant_reads.pl";
my $script_split      = "$Bin/run_split_reads.pl";
my $script_merge      = "$Bin/merge_evidence.pl";

check_script($script_depth);
check_script($script_discordant);
check_script($script_split);
check_script($script_merge);

# -----------------------------
# Prepare sample list
# -----------------------------
my @samples;

if ($bam) {
    $bam = abs_path($bam);
    push @samples, {
        sample => $sample,
        bam    => $bam,
    };
}

if ($bam_list) {
    $bam_list = abs_path($bam_list);
    die "[ERROR] BAM list file not found: $bam_list\n" unless -s $bam_list;
    @samples = read_bam_list($bam_list);
}

die "[ERROR] No valid sample found\n" unless @samples;

# -----------------------------
# Print pipeline information
# -----------------------------
print "\n";
print "============================================================\n";
print " HCMExonDel pipeline started\n";
print "================================================------------\n";
print "Config file : $config\n";
print "Target BED  : $target_bed\n";
print "Outdir      : $outdir\n";
print "Threads     : $threads\n";
print "Samples     : " . scalar(@samples) . "\n";
print "============================================================\n\n";

# -----------------------------
# Run pipeline for each sample
# -----------------------------
foreach my $item (@samples) {
    my $sid = $item->{sample};
    my $in_bam = $item->{bam};

    print "\n";
    print "------------------------------------------------------------\n";
    print "[INFO] Processing sample: $sid\n";
    print "------------------------------------------------------------\n";

    run_one_sample(
        sample     => $sid,
        bam        => $in_bam,
        config     => $config,
        target_bed => $target_bed,
        outdir     => $outdir,
        threads    => $threads,
        keep_tmp   => $keep_tmp,
        force      => $force,
    );
}

print "\n";
print "============================================================\n";
print " HCMExonDel pipeline finished successfully\n";
print "============================================================\n\n";

exit 0;


# ============================================================
# Subroutines
# ============================================================

sub run_one_sample {
    my %args = @_;

    my $sid        = $args{sample};
    my $bam        = $args{bam};
    my $config     = $args{config};
    my $target_bed = $args{target_bed};
    my $base_out   = $args{outdir};
    my $threads    = $args{threads};
    my $keep_tmp   = $args{keep_tmp};
    my $force      = $args{force};

    check_bam($bam);

    my $sample_outdir = "$base_out/$sid";

    my $log_dir        = "$sample_outdir/00.log";
    my $depth_dir      = "$sample_outdir/01.depth";
    my $discordant_dir = "$sample_outdir/02.discordant_reads";
    my $split_dir      = "$sample_outdir/03.split_reads";
    my $candidate_dir  = "$sample_outdir/04.candidates";
    my $report_dir     = "$sample_outdir/05.report";
    my $tmp_dir        = "$sample_outdir/tmp";

    if (-d $sample_outdir && !$force) {
        die "[ERROR] Output directory already exists: $sample_outdir\n"
          . "        Use --force to overwrite or choose another --outdir\n";
    }

    make_path($log_dir);
    make_path($depth_dir);
    make_path($discordant_dir);
    make_path($split_dir);
    make_path($candidate_dir);
    make_path($report_dir);
    make_path($tmp_dir);

    my $pipeline_log = "$log_dir/$sid.pipeline.log";

    log_msg($pipeline_log, "Sample: $sid");
    log_msg($pipeline_log, "BAM: $bam");
    log_msg($pipeline_log, "Output directory: $sample_outdir");

    my $depth_out      = "$depth_dir/$sid.depth_candidates.tsv";
    my $discordant_out = "$discordant_dir/$sid.discordant_reads.tsv";
    my $split_out      = "$split_dir/$sid.split_reads.tsv";
    my $final_out      = "$report_dir/$sid.final_report.tsv";

    # -----------------------------
    # Step 1: gene/exon depth
    # -----------------------------
    my $cmd_depth = join(" ",
        "perl", shell_quote($script_depth),
        "--config", shell_quote($config),
        "--bam", shell_quote($bam),
        "--sample", shell_quote($sid),
        "--target-bed", shell_quote($target_bed),
        "--out", shell_quote($depth_out),
        "--threads", $threads
    );

    run_cmd($cmd_depth, $pipeline_log, "Step 1: gene/exon depth analysis");

    # -----------------------------
    # Step 2: discordant reads
    # -----------------------------
    my $cmd_discordant = join(" ",
        "perl", shell_quote($script_discordant),
        "--config", shell_quote($config),
        "--bam", shell_quote($bam),
        "--sample", shell_quote($sid),
        "--candidate", shell_quote($depth_out),
        "--out", shell_quote($discordant_out),
        "--threads", $threads
    );

    run_cmd($cmd_discordant, $pipeline_log, "Step 2: discordant reads analysis");

    # -----------------------------
    # Step 3: split reads
    # -----------------------------
    my $cmd_split = join(" ",
        "perl", shell_quote($script_split),
        "--config", shell_quote($config),
        "--bam", shell_quote($bam),
        "--sample", shell_quote($sid),
        "--candidate", shell_quote($depth_out),
        "--out", shell_quote($split_out),
        "--threads", $threads
    );

    run_cmd($cmd_split, $pipeline_log, "Step 3: split reads analysis");

    # -----------------------------
    # Step 4: merge evidence
    # -----------------------------
    my $cmd_merge = join(" ",
        "perl", shell_quote($script_merge),
        "--config", shell_quote($config),
        "--sample", shell_quote($sid),
        "--depth", shell_quote($depth_out),
        "--discordant", shell_quote($discordant_out),
        "--split", shell_quote($split_out),
        "--out", shell_quote($final_out)
    );

    run_cmd($cmd_merge, $pipeline_log, "Step 4: merge evidence");

    # -----------------------------
    # Clean tmp directory
    # -----------------------------
    if (!$keep_tmp && -d $tmp_dir) {
        system("rm -rf " . shell_quote($tmp_dir));
        log_msg($pipeline_log, "Temporary directory removed: $tmp_dir");
    }

    log_msg($pipeline_log, "Sample finished: $sid");
    log_msg($pipeline_log, "Final report: $final_out");

    print "[INFO] Finished sample: $sid\n";
    print "[INFO] Final report: $final_out\n";
}


sub read_config {
    my ($file) = @_;

    my %conf;

    open my $fh, "<", $file or die "[ERROR] Cannot open config file: $file\n";

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        if ($line =~ /^\s*([^=\s]+)\s*=\s*(.*?)\s*$/) {
            my $key = $1;
            my $val = $2;

            $val =~ s/^['"]//;
            $val =~ s/['"]$//;

            $conf{$key} = $val;
        }
    }

    close $fh;

    return %conf;
}


sub read_bam_list {
    my ($file) = @_;

    my @list;

    open my $fh, "<", $file or die "[ERROR] Cannot open BAM list: $file\n";

    my $line_no = 0;

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        $line_no++;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        my @fields = split /\t/, $line;

        # Skip header line
        if ($line_no == 1 && $fields[0] =~ /^SampleID$/i) {
            next;
        }

        die "[ERROR] Invalid BAM list format at line $line_no: $line\n"
            unless @fields >= 2;

        my $sid = $fields[0];
        my $bam = $fields[1];

        die "[ERROR] Empty SampleID at line $line_no\n" unless $sid;
        die "[ERROR] Empty BAM path at line $line_no\n" unless $bam;

        $bam = abs_path($bam);
        die "[ERROR] BAM file not found at line $line_no: $fields[1]\n"
            unless $bam && -s $bam;

        push @list, {
            sample => $sid,
            bam    => $bam,
        };
    }

    close $fh;

    return @list;
}


sub check_bam {
    my ($bam) = @_;

    die "[ERROR] BAM file not found: $bam\n" unless -s $bam;

    my $bai1 = "$bam.bai";
    my $bai2 = $bam;
    $bai2 =~ s/\.bam$/.bai/;

    die "[ERROR] BAM index not found for: $bam\n"
        unless (-s $bai1 || -s $bai2);

    return 1;
}


sub check_script {
    my ($script) = @_;

    die "[ERROR] Required script not found: $script\n" unless -s $script;

    return 1;
}


sub run_cmd {
    my ($cmd, $log, $step_name) = @_;

    print "[INFO] $step_name\n";
    print "[CMD]  $cmd\n";

    log_msg($log, "------------------------------------------------------------");
    log_msg($log, $step_name);
    log_msg($log, "CMD: $cmd");

    my $start_time = scalar localtime();
    log_msg($log, "Start time: $start_time");

    my $ret = system($cmd);

    if ($ret != 0) {
        my $exit_code = $ret >> 8;
        log_msg($log, "[ERROR] Command failed with exit code: $exit_code");
        die "[ERROR] $step_name failed. Please check log: $log\n";
    }

    my $end_time = scalar localtime();
    log_msg($log, "End time: $end_time");
    log_msg($log, "Status: success");
}


sub log_msg {
    my ($log, $msg) = @_;

    open my $fh, ">>", $log or die "[ERROR] Cannot write log file: $log\n";
    print $fh "[" . scalar(localtime()) . "] $msg\n";
    close $fh;
}


sub shell_quote {
    my ($str) = @_;

    return "''" unless defined $str;

    $str =~ s/'/'"'"'/g;
    return "'$str'";
}


sub usage {
    return <<"USAGE";

Usage:

  Single-sample mode:

    perl bin/hcm_exondel.pl \\
      --config conf/hcm_exondel.conf \\
      --bam sample.sorted.bam \\
      --sample SAMPLE001 \\
      --outdir results

  Batch mode:

    perl bin/hcm_exondel.pl \\
      --config conf/hcm_exondel.conf \\
      --bam-list example/input_bam.list \\
      --outdir results

Required arguments:

  --config       Config file
  --bam          Coordinate-sorted BAM file
  --sample       Sample ID, required when using --bam
  --bam-list     BAM list file, tab-delimited, with columns: SampleID BAM

Optional arguments:

  --outdir       Output directory [default: results]
  --threads      Number of threads [default: config THREADS or 4]
  --target-bed   Target exon BED file
  --force        Overwrite existing sample output directory
  --keep-tmp     Keep temporary files
  --help         Show this help message

Input BAM requirement:

  1. BAM must be coordinate-sorted.
  2. BAM index must exist.
  3. BAM chromosome names should be consistent with target BED.
  4. BAM should be generated from WGS data.

Output structure:

  results/
  └── SAMPLE001/
      ├── 00.log/
      ├── 01.depth/
      ├── 02.discordant_reads/
      ├── 03.split_reads/
      ├── 04.candidates/
      └── 05.report/

Core scripts called by this pipeline:

  bin/run_gene_mean_depth.pl
  bin/run_discordant_reads.pl
  bin/run_split_reads.pl
  bin/merge_evidence.pl

USAGE
}

