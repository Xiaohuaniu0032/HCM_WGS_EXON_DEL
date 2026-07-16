#!/usr/bin/env perl
# -*- coding: utf-8 -*-
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path);

# ============================================================
# run_sample_base_depth_ratio.pl
#
# Function:
#   Sample-level wrapper for base-level depth ratio calculation.
#
# Input:
#   SAMPLE.depth_candidates.gene_depth_files/*.depth.tsv
#
# For each per-gene depth file:
#   call calc_depth_ratio.pl
#
# Per-gene outputs:
#   SAMPLE.depth_candidates.base_depth_ratio/*.base_depth_ratio.tsv
#   SAMPLE.depth_candidates.base_depth_ratio/*.base_depth.summary.tsv
#
# Optional merged outputs:
#   SAMPLE.depth_candidates.base_depth_ratio.all.tsv
#   SAMPLE.depth_candidates.base_depth.summary.all.tsv
#
# Notes:
#   1. This script does not read config directly.
#   2. HCMExonDel.pl controls whether this step is run.
#   3. HCMExonDel.pl also controls whether merged all.tsv files
#      are generated.
#   4. By default, this script only generates per-gene files.
#   5. The merged ratio all.tsv can be very large for WGS data.
# ============================================================

my %opt;

GetOptions(
    "in-dir=s"      => \$opt{in_dir},
    "out-dir=s"     => \$opt{out_dir},
    "perl=s"        => \$opt{perl},
    "script=s"      => \$opt{script},
    "ratio-all=s"   => \$opt{ratio_all},
    "summary-all=s" => \$opt{summary_all},
    "help|h"        => \$opt{help},
) or die usage();

if ($opt{help}) {
    print usage();
    exit 0;
}

for my $k (qw/in_dir out_dir perl script/) {
    die "[ERROR] Missing required option --" . option_name($k) . "\n" . usage()
        unless defined $opt{$k} && $opt{$k} ne "";
}

if ((defined $opt{ratio_all} && $opt{ratio_all} ne "") xor
    (defined $opt{summary_all} && $opt{summary_all} ne "")) {
    die "[ERROR] --ratio-all and --summary-all must be provided together.\n";
}

die "[ERROR] Input gene depth directory not found: $opt{in_dir}\n"
    unless -d $opt{in_dir};

die "[ERROR] calc_depth_ratio.pl not found: $opt{script}\n"
    unless -s $opt{script};

check_executable($opt{perl}, "perl");

$opt{in_dir} = abs_path($opt{in_dir});
$opt{script} = abs_path($opt{script});

make_path($opt{out_dir}) unless -d $opt{out_dir};
$opt{out_dir} = abs_path($opt{out_dir});

if (defined $opt{ratio_all} && $opt{ratio_all} ne "") {
    my $ratio_parent = dirname($opt{ratio_all});
    make_path($ratio_parent) if defined $ratio_parent && $ratio_parent ne "." && !-d $ratio_parent;

    my $summary_parent = dirname($opt{summary_all});
    make_path($summary_parent) if defined $summary_parent && $summary_parent ne "." && !-d $summary_parent;
}

log_msg("Sample-level base depth ratio calculation started");
log_msg("Input depth dir       : $opt{in_dir}");
log_msg("Output ratio dir      : $opt{out_dir}");
log_msg("Perl                  : $opt{perl}");
log_msg("calc_depth_ratio.pl   : $opt{script}");

if (defined $opt{ratio_all} && $opt{ratio_all} ne "") {
    log_msg("Merged ratio output   : $opt{ratio_all}");
    log_msg("Merged summary output : $opt{summary_all}");
}
else {
    log_msg("Merged all outputs    : disabled");
}

my @depth_files = sort glob("$opt{in_dir}/*.depth.tsv");

die "[ERROR] No *.depth.tsv files found in $opt{in_dir}\n"
    unless @depth_files;

log_msg("Input depth files found: " . scalar(@depth_files));

my @ratio_files;
my @summary_files;

for my $depth_file (@depth_files) {
    my $base = basename($depth_file);

    my $ratio_base = $base;
    $ratio_base =~ s/\.depth\.tsv$/.base_depth_ratio.tsv/;
    $ratio_base .= ".base_depth_ratio.tsv" if $ratio_base eq $base;

    my $summary_base = $base;
    $summary_base =~ s/\.depth\.tsv$/.base_depth.summary.tsv/;
    $summary_base .= ".base_depth.summary.tsv" if $summary_base eq $base;

    my $ratio_out   = "$opt{out_dir}/$ratio_base";
    my $summary_out = "$opt{out_dir}/$summary_base";

    log_msg("Calculating base depth ratio: $base");

    my @cmd = (
        $opt{perl},
        $opt{script},
        $depth_file,
        $ratio_out,
        $summary_out,
    );

    run_cmd(@cmd);

    die "[ERROR] Ratio output not generated: $ratio_out\n"
        unless -s $ratio_out;

    die "[ERROR] Summary output not generated: $summary_out\n"
        unless -s $summary_out;

    push @ratio_files,   $ratio_out;
    push @summary_files, $summary_out;
}

if (defined $opt{ratio_all} && $opt{ratio_all} ne "") {
    merge_ratio_files(\@ratio_files, $opt{ratio_all});
    merge_summary_files(\@summary_files, $opt{summary_all});

    log_msg("Merged ratio output written   : $opt{ratio_all}");
    log_msg("Merged summary output written : $opt{summary_all}");
}

log_msg("Per-gene ratio files generated  : " . scalar(@ratio_files));
log_msg("Per-gene summary files generated: " . scalar(@summary_files));
log_msg("Sample-level base depth ratio calculation finished");

exit 0;

# ============================================================
# Merge per-gene ratio outputs
# ============================================================

sub merge_ratio_files {
    my ($files_ref, $out_file) = @_;

    open my $OUT, ">", $out_file
        or die "[ERROR] Cannot write merged ratio output $out_file: $!\n";

    my $header_written = 0;
    my $total_records  = 0;

    for my $file (@$files_ref) {
        open my $IN, "<", $file
            or die "[ERROR] Cannot open ratio file $file: $!\n";

        my $header = <$IN>;

        if (!defined $header) {
            close $IN;
            next;
        }

        chomp $header;
        $header =~ s/\r$//;

        my @cols = split /\t/, $header, -1;

        my $expected = join(
            "\t",
            "#CHROM",
            qw/
              POS
              Depth
              Median_Depth
              Mean_Depth
              Ratio_By_Median
              Ratio_By_Mean
            /
        );

        my $observed = join("\t", @cols);

        die "[ERROR] Unexpected ratio header in $file\n"
          . "[ERROR] Expected: $expected\n"
          . "[ERROR] Observed: $observed\n"
            unless $observed eq $expected;

        if (!$header_written) {
            print $OUT join("\t", "Depth_File", @cols), "\n";
            $header_written = 1;
        }

        my $depth_file_name = basename($file);
        $depth_file_name =~ s/\.base_depth_ratio\.tsv$/.depth.tsv/;

        while (my $line = <$IN>) {
            chomp $line;
            $line =~ s/\r$//;

            next if $line =~ /^\s*$/;

            print $OUT $depth_file_name, "\t", $line, "\n";
            $total_records++;
        }

        close $IN;
    }

    close $OUT;

    log_msg("Merged base depth ratio records: $total_records");
}

sub merge_summary_files {
    my ($files_ref, $out_file) = @_;

    open my $OUT, ">", $out_file
        or die "[ERROR] Cannot write merged summary output $out_file: $!\n";

    print $OUT join("\t", qw/Depth_File Metric Value/), "\n";

    my $total_records = 0;

    for my $file (@$files_ref) {
        open my $IN, "<", $file
            or die "[ERROR] Cannot open summary file $file: $!\n";

        my $header = <$IN>;

        if (!defined $header) {
            close $IN;
            next;
        }

        chomp $header;
        $header =~ s/\r$//;

        my @cols = split /\t/, $header, -1;
        my $observed = join("\t", @cols);
        my $expected = "Metric\tValue";

        die "[ERROR] Unexpected summary header in $file\n"
          . "[ERROR] Expected: $expected\n"
          . "[ERROR] Observed: $observed\n"
            unless $observed eq $expected;

        my $depth_file_name = basename($file);
        $depth_file_name =~ s/\.base_depth\.summary\.tsv$/.depth.tsv/;

        while (my $line = <$IN>) {
            chomp $line;
            $line =~ s/\r$//;

            next if $line =~ /^\s*$/;

            print $OUT $depth_file_name, "\t", $line, "\n";
            $total_records++;
        }

        close $IN;
    }

    close $OUT;

    log_msg("Merged base depth summary records: $total_records");
}

# ============================================================
# General helpers
# ============================================================

sub check_executable {
    my ($cmd, $name) = @_;

    if ($cmd =~ /\//) {
        die "[ERROR] $name executable not found or not executable: $cmd\n"
            unless -x $cmd;
    }
    else {
        my $ret = system("command -v $cmd >/dev/null 2>&1");
        die "[ERROR] $name executable not found in PATH: $cmd\n"
            if $ret != 0;
    }

    return 1;
}

sub run_cmd {
    my (@cmd) = @_;

    print STDERR "[CMD] ", join(" ", map { shell_quote($_) } @cmd), "\n";

    system(@cmd) == 0
        or die "[ERROR] Command failed: " . join(" ", @cmd) . "\n";
}

sub shell_quote {
    my ($s) = @_;

    return "''" if !defined $s || $s eq "";

    if ($s =~ /^[A-Za-z0-9_\.\-\/:=]+$/) {
        return $s;
    }

    $s =~ s/'/'\\''/g;

    return "'$s'";
}

sub option_name {
    my ($k) = @_;

    $k =~ s/_/-/g;

    return $k;
}

sub log_msg {
    my ($msg) = @_;

    print STDERR "[INFO] $msg\n";
}

sub usage {
    return <<'USAGE';
Usage:
  perl run_sample_base_depth_ratio.pl \
    --in-dir SAMPLE.depth_candidates.gene_depth_files \
    --out-dir SAMPLE.depth_candidates.base_depth_ratio \
    --perl /usr/bin/perl \
    --script bin/calc_depth_ratio.pl

Required:
  --in-dir
      Directory containing per-gene depth files:
        *.depth.tsv

  --out-dir
      Output directory for per-gene base depth ratio files:
        *.base_depth_ratio.tsv
        *.base_depth.summary.tsv

  --perl
      Perl executable passed from HCMExonDel.pl.

  --script
      Path to bin/calc_depth_ratio.pl.

Optional:
  --ratio-all
      Merged sample-level base depth ratio file:
        SAMPLE.depth_candidates.base_depth_ratio.all.tsv

  --summary-all
      Merged sample-level base depth summary file:
        SAMPLE.depth_candidates.base_depth.summary.all.tsv

      Note:
        --ratio-all and --summary-all must be provided together.

Output:
  1. Per-gene ratio files:
       out-dir/*.base_depth_ratio.tsv

  2. Per-gene summary files:
       out-dir/*.base_depth.summary.tsv

  3. Optional merged ratio file:
       --ratio-all

  4. Optional merged summary file:
       --summary-all

Example without merged all.tsv:
  perl bin/run_sample_base_depth_ratio.pl \
    --in-dir test/test_results/25B09089386/01.depth/25B09089386.depth_candidates.gene_depth_files \
    --out-dir test/test_results/25B09089386/01.depth/25B09089386.depth_candidates.base_depth_ratio \
    --perl /usr/bin/perl \
    --script bin/calc_depth_ratio.pl

Example with merged all.tsv:
  perl bin/run_sample_base_depth_ratio.pl \
    --in-dir test/test_results/25B09089386/01.depth/25B09089386.depth_candidates.gene_depth_files \
    --out-dir test/test_results/25B09089386/01.depth/25B09089386.depth_candidates.base_depth_ratio \
    --perl /usr/bin/perl \
    --script bin/calc_depth_ratio.pl \
    --ratio-all test/test_results/25B09089386/01.depth/25B09089386.depth_candidates.base_depth_ratio.all.tsv \
    --summary-all test/test_results/25B09089386/01.depth/25B09089386.depth_candidates.base_depth.summary.all.tsv
USAGE
}

