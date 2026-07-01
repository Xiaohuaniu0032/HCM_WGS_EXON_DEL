#!/usr/bin/env perl
# -*- coding: utf-8 -*-
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path);

# ============================================================
# run_sample_cusum_depth.pl
#
# Function:
#   Sample-level wrapper for CUSUM depth deletion calling.
#
# Input:
#   SAMPLE.depth_candidates.gene_depth_files/*.depth.tsv
#
# For each per-gene depth file:
#   call cusum_depth_del.pl
#
# Output:
#   1. Per-gene CUSUM files:
#        SAMPLE.depth_candidates.cusum/*.cusum.tsv
#
#   2. Merged sample-level CUSUM file:
#        SAMPLE.depth_candidates.cusum.all.tsv
#
# Output format of merged file:
#   Depth_File
#   CHROM
#   START
#   END
#   LENGTH
#   N_BINS
#   BASELINE
#   MEAN_DEPTH
#   MEAN_RATIO
#   MIN_RATIO
#   MAX_CUSUM
#   TYPE
#
# Notes:
#   1. This script does not define algorithm defaults.
#   2. All CUSUM parameters must be explicitly provided in config.
#   3. CUSUM_BASELINE can be "auto" or a numeric value.
#   4. Coordinates are inherited from cusum_depth_del.pl.
# ============================================================

my %opt;

GetOptions(
    "config|conf=s" => \$opt{config},
    "in-dir=s"      => \$opt{in_dir},
    "out-dir=s"     => \$opt{out_dir},
    "out=s"         => \$opt{out},
    "perl=s"        => \$opt{perl},
    "script=s"      => \$opt{script},
    "help|h"        => \$opt{help},
) or die usage();

if ($opt{help}) {
    print usage();
    exit 0;
}

for my $k (qw/config in_dir out_dir out perl script/) {
    die "[ERROR] Missing required option --" . option_name($k) . "\n" . usage()
        unless defined $opt{$k} && $opt{$k} ne "";
}

die "[ERROR] Config file not found: $opt{config}\n"
    unless -s $opt{config};

die "[ERROR] Input gene depth directory not found: $opt{in_dir}\n"
    unless -d $opt{in_dir};

die "[ERROR] CUSUM script not found: $opt{script}\n"
    unless -s $opt{script};

check_executable($opt{perl}, "perl");

$opt{config} = abs_path($opt{config});
$opt{in_dir} = abs_path($opt{in_dir});
$opt{script} = abs_path($opt{script});

make_path($opt{out_dir}) unless -d $opt{out_dir};
$opt{out_dir} = abs_path($opt{out_dir});

my $out_parent = dirname($opt{out});
make_path($out_parent) if defined $out_parent && $out_parent ne "." && !-d $out_parent;

my %conf = read_config($opt{config});

my $baseline     = get_conf_required(\%conf, "CUSUM_BASELINE");
my $bin_size     = get_conf_required(\%conf, "CUSUM_BIN_SIZE");
my $k            = get_conf_required(\%conf, "CUSUM_K");
my $h            = get_conf_required(\%conf, "CUSUM_H");
my $del_ratio    = get_conf_required(\%conf, "CUSUM_DEL_RATIO");
my $min_bins     = get_conf_required(\%conf, "CUSUM_MIN_BINS");
my $min_len      = get_conf_required(\%conf, "CUSUM_MIN_LEN");
my $edge_ratio   = get_conf_required(\%conf, "CUSUM_EDGE_RATIO");
my $recover_ratio = get_conf_required(\%conf, "CUSUM_RECOVER_RATIO");
my $recover_bins = get_conf_required(\%conf, "CUSUM_RECOVER_BINS");

validate_baseline("CUSUM_BASELINE", $baseline);
validate_positive_integer("CUSUM_BIN_SIZE", $bin_size);
validate_nonnegative_number("CUSUM_K", $k);
validate_positive_number("CUSUM_H", $h);
validate_fraction("CUSUM_DEL_RATIO", $del_ratio);
validate_positive_integer("CUSUM_MIN_BINS", $min_bins);
validate_positive_integer("CUSUM_MIN_LEN", $min_len);
validate_fraction("CUSUM_EDGE_RATIO", $edge_ratio);
validate_fraction("CUSUM_RECOVER_RATIO", $recover_ratio);
validate_positive_integer("CUSUM_RECOVER_BINS", $recover_bins);

log_msg("Sample-level CUSUM depth calling started");
log_msg("Config             : $opt{config}");
log_msg("Input depth dir    : $opt{in_dir}");
log_msg("Output CUSUM dir   : $opt{out_dir}");
log_msg("Merged output      : $opt{out}");
log_msg("Perl               : $opt{perl}");
log_msg("CUSUM script       : $opt{script}");
log_msg("CUSUM_BASELINE     : $baseline");
log_msg("CUSUM_BIN_SIZE     : $bin_size");
log_msg("CUSUM_K            : $k");
log_msg("CUSUM_H            : $h");
log_msg("CUSUM_DEL_RATIO    : $del_ratio");
log_msg("CUSUM_MIN_BINS     : $min_bins");
log_msg("CUSUM_MIN_LEN      : $min_len");
log_msg("CUSUM_EDGE_RATIO   : $edge_ratio");
log_msg("CUSUM_RECOVER_RATIO: $recover_ratio");
log_msg("CUSUM_RECOVER_BINS : $recover_bins");

my @depth_files = sort glob("$opt{in_dir}/*.depth.tsv");

die "[ERROR] No *.depth.tsv files found in $opt{in_dir}\n"
    unless @depth_files;

log_msg("Input depth files found: " . scalar(@depth_files));

my @cusum_files;

for my $depth_file (@depth_files) {
    my $base = basename($depth_file);

    my $out_base = $base;
    $out_base =~ s/\.depth\.tsv$/.cusum.tsv/;

    if ($out_base eq $base) {
        $out_base .= ".cusum.tsv";
    }

    my $cusum_out = "$opt{out_dir}/$out_base";

    log_msg("Running CUSUM: $base");

    my @cmd = (
        $opt{perl},
        $opt{script},
        "--in",            $depth_file,
        "--out",           $cusum_out,
        "--bin-size",      $bin_size,
        "--baseline",      $baseline,
        "--k",             $k,
        "--h",             $h,
        "--edge-ratio",    $edge_ratio,
        "--recover-ratio", $recover_ratio,
        "--recover-bins",  $recover_bins,
        "--del-ratio",     $del_ratio,
        "--min-bins",      $min_bins,
        "--min-len",       $min_len,
    );

    run_cmd(@cmd);

    die "[ERROR] CUSUM output not generated: $cusum_out\n"
        unless -e $cusum_out;

    push @cusum_files, $cusum_out;
}

merge_cusum_files(\@cusum_files, $opt{out});

log_msg("Per-gene CUSUM files generated: " . scalar(@cusum_files));
log_msg("Merged CUSUM output written   : $opt{out}");
log_msg("Sample-level CUSUM depth calling finished");

exit 0;

# ============================================================
# Merge per-gene CUSUM outputs
# ============================================================

sub merge_cusum_files {
    my ($files_ref, $out_file) = @_;

    open my $OUT, ">", $out_file
        or die "[ERROR] Cannot write merged CUSUM output $out_file: $!\n";

    print $OUT join(
        "\t",
        qw/
          Depth_File
          CHROM
          START
          END
          LENGTH
          N_BINS
          BASELINE
          MEAN_DEPTH
          MEAN_RATIO
          MIN_RATIO
          MAX_CUSUM
          TYPE
        /
    ), "\n";

    my $total_records = 0;

    for my $file (@$files_ref) {
        open my $IN, "<", $file
            or die "[ERROR] Cannot open per-gene CUSUM output $file: $!\n";

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
            qw/
              CHROM
              START
              END
              LENGTH
              N_BINS
              BASELINE
              MEAN_DEPTH
              MEAN_RATIO
              MIN_RATIO
              MAX_CUSUM
              TYPE
            /
        );

        my $observed = join("\t", @cols);

        die "[ERROR] Unexpected CUSUM header in $file\n"
          . "[ERROR] Expected: $expected\n"
          . "[ERROR] Observed: $observed\n"
            unless $observed eq $expected;

        my $depth_file_name = basename($file);
        $depth_file_name =~ s/\.cusum\.tsv$/.depth.tsv/;

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

    log_msg("Merged CUSUM records: $total_records");
}

# ============================================================
# Config and validation
# ============================================================

sub read_config {
    my ($file) = @_;

    open my $fh, "<", $file
        or die "[ERROR] Cannot open config file $file: $!\n";

    my %conf;

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        $line =~ s/\s+#.*$//;

        next unless $line =~ /^\s*([^=\s]+)\s*=\s*(.*?)\s*$/;

        my $key = $1;
        my $val = $2;

        $val =~ s/^\s+|\s+$//g;
        $val =~ s/^["']//;
        $val =~ s/["']$//;

        $conf{$key} = $val;
    }

    close $fh;

    return %conf;
}

sub get_conf_required {
    my ($conf_ref, $key) = @_;

    die "[ERROR] Required config key '$key' is missing\n"
        unless exists $conf_ref->{$key}
        && defined $conf_ref->{$key}
        && $conf_ref->{$key} ne "";

    return $conf_ref->{$key};
}

sub validate_baseline {
    my ($name, $v) = @_;

    return 1 if defined $v && $v eq "auto";

    validate_positive_number($name, $v);

    return 1;
}

sub validate_fraction {
    my ($name, $v) = @_;

    die "[ERROR] $name must be numeric: $v\n"
        unless defined $v
        && $v =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;

    die "[ERROR] $name must be > 0 and <= 1: $v\n"
        unless $v > 0 && $v <= 1;

    return 1;
}

sub validate_positive_integer {
    my ($name, $v) = @_;

    die "[ERROR] $name must be a positive integer: $v\n"
        unless defined $v && $v =~ /^\d+$/ && $v >= 1;

    return 1;
}

sub validate_positive_number {
    my ($name, $v) = @_;

    die "[ERROR] $name must be numeric: $v\n"
        unless defined $v
        && $v =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;

    die "[ERROR] $name must be > 0: $v\n"
        unless $v > 0;

    return 1;
}

sub validate_nonnegative_number {
    my ($name, $v) = @_;

    die "[ERROR] $name must be numeric: $v\n"
        unless defined $v
        && $v =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;

    die "[ERROR] $name must be >= 0: $v\n"
        unless $v >= 0;

    return 1;
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
  perl run_sample_cusum_depth.pl \
    --config conf/hcm_exondel.example.conf \
    --in-dir SAMPLE.depth_candidates.gene_depth_files \
    --out-dir SAMPLE.depth_candidates.cusum \
    --out SAMPLE.depth_candidates.cusum.all.tsv \
    --perl /usr/bin/perl \
    --script bin/cusum_depth_del.pl

Required:
  --config
      HCMExonDel config file.

  --in-dir
      Directory containing per-gene depth files:
        *.depth.tsv

  --out-dir
      Output directory for per-gene CUSUM result files:
        *.cusum.tsv

  --out
      Merged sample-level CUSUM output:
        SAMPLE.depth_candidates.cusum.all.tsv

  --perl
      Perl executable passed from HCMExonDel.pl.

  --script
      Path to bin/cusum_depth_del.pl.

Required config keys:
  CUSUM_BASELINE
  CUSUM_BIN_SIZE
  CUSUM_K
  CUSUM_H
  CUSUM_DEL_RATIO
  CUSUM_MIN_BINS
  CUSUM_MIN_LEN
  CUSUM_EDGE_RATIO
  CUSUM_RECOVER_RATIO
  CUSUM_RECOVER_BINS

Output:
  1. Per-gene CUSUM files:
       out-dir/*.cusum.tsv

  2. Merged sample-level CUSUM file:
       --out

Merged output columns:
  Depth_File
  CHROM
  START
  END
  LENGTH
  N_BINS
  BASELINE
  MEAN_DEPTH
  MEAN_RATIO
  MIN_RATIO
  MAX_CUSUM
  TYPE

Example:
  perl bin/run_sample_cusum_depth.pl \
    --config conf/hcm_exondel.example.conf \
    --in-dir test/test_results/25B09089386/01.depth/25B09089386.depth_candidates.gene_depth_files \
    --out-dir test/test_results/25B09089386/01.depth/25B09089386.depth_candidates.cusum \
    --out test/test_results/25B09089386/01.depth/25B09089386.depth_candidates.cusum.all.tsv \
    --perl /usr/bin/perl \
    --script bin/cusum_depth_del.pl
USAGE
}

