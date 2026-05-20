#!/usr/bin/env perl
# -*- coding: utf-8 -*-

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path);

# ============================================================
# extract_candidate_gene_bam.pl
#
# Function:
#   Extract BAM files for genes listed in merged_candidates.tsv.
#
# Input:
#   --candidate merged_candidates.tsv
#     Must contain column: Gene
#
# Behavior:
#   1. If candidate file is empty, skip and exit 0.
#   2. If candidate file has only header but no valid Gene value, skip and exit 0.
#   3. If Gene column is missing, report ERROR and exit non-zero.
#   4. If valid genes exist but none match gene TXT, report ERROR and exit non-zero.
#   5. Default flank is 0, meaning strict gene-coordinate extraction.
# ============================================================

my ($config, $candidate_tsv, $bam, $outdir, $flank, $help);

$flank = 0;

GetOptions(
    "config|conf=s" => \$config,
    "candidate|c=s" => \$candidate_tsv,
    "bam|b=s"       => \$bam,
    "outdir|o=s"    => \$outdir,
    "flank=i"       => \$flank,
    "help|h"        => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die "[ERROR] --config is required\n"    unless defined $config && $config ne "";
die "[ERROR] --candidate is required\n" unless defined $candidate_tsv && $candidate_tsv ne "";
die "[ERROR] --bam is required\n"       unless defined $bam && $bam ne "";
die "[ERROR] --outdir is required\n"    unless defined $outdir && $outdir ne "";

die "[ERROR] Config file not found: $config\n" unless -s $config;
die "[ERROR] BAM file not found: $bam\n"       unless -s $bam;

die "[ERROR] --flank must be a non-negative integer\n"
    unless defined $flank && $flank =~ /^\d+$/;

$config = abs_path($config);
$bam    = abs_path($bam);

make_path($outdir) unless -d $outdir;
$outdir = abs_path($outdir);

my $project_root = detect_project_root($config);
my %conf = read_config($config);

my $gene_txt = get_required_conf_path(
    \%conf,
    "REFSEQ_MANE_SELECT_GENE_TXT",
    $project_root,
    "gene coordinate file"
);

my $samtools = get_samtools_from_conf(\%conf);
check_executable($samtools, "samtools");

print STDERR "[INFO] Extract candidate gene BAM started\n";
print STDERR "[INFO] Config        : $config\n";
print STDERR "[INFO] Project root  : $project_root\n";
print STDERR "[INFO] Candidate TSV : $candidate_tsv\n";
print STDERR "[INFO] Gene TXT      : $gene_txt\n";
print STDERR "[INFO] Input BAM     : $bam\n";
print STDERR "[INFO] Output dir    : $outdir\n";
print STDERR "[INFO] Samtools      : $samtools\n";
print STDERR "[INFO] Flank         : $flank\n";

# ------------------------------------------------------------
# Step 1. Candidate file pre-check
# ------------------------------------------------------------

if (!-e $candidate_tsv) {
    die "[ERROR] Candidate file not found: $candidate_tsv\n";
}

if (!-s $candidate_tsv) {
    print STDERR "[INFO] Candidate file is empty. Skip candidate gene BAM extraction.\n";
    exit 0;
}

$candidate_tsv = abs_path($candidate_tsv);

# ------------------------------------------------------------
# Step 2. Read Gene column from merged_candidates.tsv
# ------------------------------------------------------------

my %target_genes = read_candidate_genes($candidate_tsv);
my $gene_count = scalar keys %target_genes;

if ($gene_count == 0) {
    print STDERR "[INFO] No valid Gene value found in candidate file. Skip candidate gene BAM extraction.\n";
    exit 0;
}

print STDERR "[INFO] Candidate genes found: $gene_count\n";

# ------------------------------------------------------------
# Step 3. Read gene coordinates from RefSeq_MANE_Select.gene.txt
# ------------------------------------------------------------

my %gene_regions = read_gene_regions($gene_txt, \%target_genes, $flank);
my $matched_count = scalar keys %gene_regions;

print STDERR "[INFO] Genes matched in gene TXT: $matched_count / $gene_count\n";

for my $gene (sort keys %target_genes) {
    if (!exists $gene_regions{$gene}) {
        print STDERR "[WARN] Gene not found in gene TXT: $gene\n";
    }
}

die "[ERROR] No candidate gene matched coordinates in $gene_txt\n"
    if $matched_count == 0;

# ------------------------------------------------------------
# Step 4. Extract BAM and generate index
# ------------------------------------------------------------

for my $gene (sort keys %gene_regions) {
    my @regions = @{ $gene_regions{$gene} };

    my $safe_gene = sanitize_filename($gene);
    my $out_bam   = "$outdir/$safe_gene.bam";

    print STDERR "[INFO] Extracting gene: $gene\n";
    print STDERR "[INFO] Regions       : ", join(",", @regions), "\n";
    print STDERR "[INFO] Output BAM    : $out_bam\n";

    my @view_cmd = (
        $samtools,
        "view",
        "-b",
        "-o",
        $out_bam,
        $bam,
        @regions
    );

    run_cmd(@view_cmd);

    die "[ERROR] Output BAM not generated or empty: $out_bam\n"
        unless -s $out_bam;

    my @index_cmd = (
        $samtools,
        "index",
        $out_bam
    );

    run_cmd(@index_cmd);

    die "[ERROR] BAM index not generated: $out_bam.bai\n"
        unless -s "$out_bam.bai";

    print STDERR "[INFO] Finished gene: $gene\n";
}

print STDERR "[INFO] All done.\n";

exit 0;

# ============================================================
# Subroutines
# ============================================================

sub read_config {
    my ($file) = @_;

    open my $fh, "<", $file or die "[ERROR] Cannot open config file $file: $!\n";

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

sub detect_project_root {
    my ($config_path) = @_;

    my $dir = dirname($config_path);

    if ($dir =~ /\/conf$/) {
        my $root = dirname($dir);
        return abs_path($root);
    }

    return abs_path($dir);
}

sub get_required_conf_path {
    my ($conf_ref, $key, $project_root, $desc) = @_;

    die "[ERROR] Required config key '$key' not found for $desc\n"
        unless exists $conf_ref->{$key}
            && defined $conf_ref->{$key}
            && $conf_ref->{$key} ne "";

    my $path = $conf_ref->{$key};

    if ($path !~ /^\//) {
        $path = "$project_root/$path";
    }

    $path = abs_path($path) if -e $path;

    die "[ERROR] $desc not found: $path\n"
        unless defined $path && -s $path;

    return $path;
}

sub get_samtools_from_conf {
    my ($conf_ref) = @_;

    my $samtools = "samtools";

    if (exists $conf_ref->{SAMTOOLS} && $conf_ref->{SAMTOOLS} ne "") {
        $samtools = $conf_ref->{SAMTOOLS};
    }
    elsif (exists $conf_ref->{SAMTOOLS_PATH} && $conf_ref->{SAMTOOLS_PATH} ne "") {
        $samtools = $conf_ref->{SAMTOOLS_PATH};
    }

    return $samtools;
}

sub check_executable {
    my ($cmd, $name) = @_;

    if ($cmd =~ /\//) {
        die "[ERROR] $name not executable: $cmd\n" unless -x $cmd;
    }
    else {
        my $check = system("command -v $cmd >/dev/null 2>&1");
        die "[ERROR] $name not found in PATH: $cmd\n" if $check != 0;
    }
}

sub read_candidate_genes {
    my ($file) = @_;

    open my $fh, "<", $file or die "[ERROR] Cannot open candidate file $file: $!\n";

    my $header = <$fh>;

    if (!defined $header) {
        close $fh;
        return ();
    }

    chomp $header;
    $header =~ s/\r$//;

    my @cols = split /\t/, $header, -1;

    my $gene_idx = -1;

    for my $i (0 .. $#cols) {
        if ($cols[$i] eq "Gene") {
            $gene_idx = $i;
            last;
        }
    }

    die "[ERROR] Required column 'Gene' not found in candidate file: $file\n"
        if $gene_idx < 0;

    my %genes;

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;

        my @f = split /\t/, $line, -1;

        next unless defined $f[$gene_idx];

        my $gene_field = $f[$gene_idx];

        $gene_field =~ s/^\s+|\s+$//g;

        next if $gene_field eq "";
        next if $gene_field eq ".";
        next if uc($gene_field) eq "NA";

        my @items = split /[;,|]/, $gene_field;

        for my $gene (@items) {
            $gene =~ s/^\s+|\s+$//g;

            next if $gene eq "";
            next if $gene eq ".";
            next if uc($gene) eq "NA";

            $genes{$gene} = 1;
        }
    }

    close $fh;

    return %genes;
}

sub read_gene_regions {
    my ($file, $target_genes_ref, $flank) = @_;

    open my $fh, "<", $file or die "[ERROR] Cannot open gene TXT file $file: $!\n";

    my $header = <$fh>;

    die "[ERROR] Empty gene TXT file: $file\n" unless defined $header;

    chomp $header;
    $header =~ s/\r$//;

    my @cols = split /\t/, $header, -1;

    my $gene_idx  = find_col(\@cols, qw(Gene GeneName gene gene_name gene_symbol Symbol symbol));
    my $chr_idx   = find_col(\@cols, qw(Chr Chrom Chromosome chr chrom chromosome));
    my $start_idx = find_col(\@cols, qw(Start Gene_start GeneStart start gene_start geneStart txStart));
    my $end_idx   = find_col(\@cols, qw(End Gene_end GeneEnd end gene_end geneEnd txEnd));

    if ($gene_idx < 0 || $chr_idx < 0 || $start_idx < 0 || $end_idx < 0) {
        die "[ERROR] Required columns not found in gene TXT: $file\n"
          . "        Required columns include Gene, Chrom, Start and End.\n";
    }

    my %gene_regions;

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;
        next if $line =~ /^#/;

        my @f = split /\t/, $line, -1;

        next unless defined $f[$gene_idx];
        next unless defined $f[$chr_idx];
        next unless defined $f[$start_idx];
        next unless defined $f[$end_idx];

        my $gene  = $f[$gene_idx];
        my $chr   = $f[$chr_idx];
        my $start = $f[$start_idx];
        my $end   = $f[$end_idx];

        $gene  =~ s/^\s+|\s+$//g;
        $chr   =~ s/^\s+|\s+$//g;
        $start =~ s/^\s+|\s+$//g;
        $end   =~ s/^\s+|\s+$//g;

        next unless exists $target_genes_ref->{$gene};

        unless ($start =~ /^\d+$/ && $end =~ /^\d+$/) {
            print STDERR "[WARN] Invalid coordinate skipped: $line\n";
            next;
        }

        my $region_start = $start - $flank;
        my $region_end   = $end + $flank;

        $region_start = 1 if $region_start < 1;

        if ($region_start > $region_end) {
            print STDERR "[WARN] Invalid region skipped for gene $gene: $chr:$region_start-$region_end\n";
            next;
        }

        my $region = "$chr:$region_start-$region_end";

        push @{ $gene_regions{$gene} }, $region;
    }

    close $fh;

    for my $gene (keys %gene_regions) {
        my %seen;
        my @uniq;

        for my $region (@{ $gene_regions{$gene} }) {
            next if $seen{$region}++;
            push @uniq, $region;
        }

        $gene_regions{$gene} = \@uniq;
    }

    return %gene_regions;
}

sub find_col {
    my ($cols_ref, @names) = @_;

    my %wanted = map { $_ => 1 } @names;

    for my $i (0 .. $#$cols_ref) {
        my $col = $cols_ref->[$i];
        $col =~ s/^\s+|\s+$//g;

        return $i if exists $wanted{$col};
    }

    return -1;
}

sub sanitize_filename {
    my ($name) = @_;

    $name =~ s/^\s+|\s+$//g;
    $name =~ s/[\/\\:\*\?"<>\|\s]+/_/g;

    return $name;
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

sub usage {
    return <<"USAGE";
Usage:
  perl extract_candidate_gene_bam.pl \\
    --config conf/hcm_exondel.example.conf \\
    --candidate merged_candidates.tsv \\
    --bam sample.bam \\
    --outdir gene_bam

Required:
  --config|-conf      Config file.
                      Required config key:
                      REFSEQ_MANE_SELECT_GENE_TXT

  --candidate|-c      merged_candidates.tsv file.
                      Must contain column: Gene

  --bam|-b            Input BAM file.

  --outdir|-o         Output directory.

Optional:
  --flank             Flanking size added to both sides of gene region.
                      Default: 0

  --help|-h           Show help message.

Behavior:
  If candidate file is empty or has no valid Gene value,
  the script will print INFO message and exit 0.

  If Gene column is missing, the script will report ERROR.

Example:
  perl bin/extract_candidate_gene_bam.pl \\
    --config conf/hcm_exondel.example.conf \\
    --candidate test/test_results/sample/04.candidates/sample.merged_candidates.tsv \\
    --bam /path/to/sample.bam \\
    --outdir test/test_results/sample/05.gene_bam

Output:
  gene_bam/
    FHOD3.bam
    FHOD3.bam.bai
    MYH7.bam
    MYH7.bam.bai
USAGE
}

