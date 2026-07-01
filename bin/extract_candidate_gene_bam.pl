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
#   Extract per-gene BAM files from a coordinate-sorted WGS BAM.
#
# Gene sources:
#   1. Genes from merged_candidates.tsv
#      - supported columns:
#          Gene
#          Annotated_Gene
#
#   2. Genes from HCM_CORE_GENE_LIST in config
#      - expected format:
#          Gene    Classification
#          PRKAG2  Definitive
#          ACTN2   Definitive
#      - only the first column is used as gene symbol.
#
# Final target genes:
#   union(candidate genes, HCM core genes)
#
# Coordinate source:
#   REFSEQ_MANE_SELECT_GENE_TXT in config
#
# Extraction rule:
#   Extract BAM according to gene coordinates.
#   Optional upstream/downstream flanking region can be added by --flank.
#
# Output:
#   outdir/
#     FHOD3.bam
#     FHOD3.bam.bai
#     MYH7.bam
#     MYH7.bam.bai
#
# Notes:
#   1. Gene extraction is based on full gene coordinates.
#   2. If --flank 0 is used, behavior is the same as the old version.
#   3. If candidate file is empty, the script still uses HCM_CORE_GENE_LIST.
#   4. If one gene has multiple coordinate records, all regions are extracted
#      into the same gene BAM.
# ============================================================

my ($config, $candidate_tsv, $bam, $outdir, $flank, $help);

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

$flank = 0 unless defined $flank;

die "[ERROR] --flank must be a non-negative integer: $flank\n"
    unless $flank =~ /^\d+$/ && $flank >= 0;

die "[ERROR] Config file not found: $config\n" unless -s $config;
die "[ERROR] BAM file not found: $bam\n"       unless -s $bam;

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
    "RefSeq MANE gene coordinate file"
);

my $core_gene_list = get_required_conf_path(
    \%conf,
    "HCM_CORE_GENE_LIST",
    $project_root,
    "HCM core gene list"
);

my $samtools = get_samtools_from_conf(\%conf);
check_executable($samtools, "samtools");

print STDERR "[INFO] Extract gene BAM started\n";
print STDERR "[INFO] Config             : $config\n";
print STDERR "[INFO] Project root       : $project_root\n";
print STDERR "[INFO] Candidate TSV      : $candidate_tsv\n";
print STDERR "[INFO] HCM core gene list : $core_gene_list\n";
print STDERR "[INFO] Gene coordinate TXT: $gene_txt\n";
print STDERR "[INFO] Input BAM          : $bam\n";
print STDERR "[INFO] Output dir         : $outdir\n";
print STDERR "[INFO] Samtools           : $samtools\n";
print STDERR "[INFO] Flank size         : $flank bp\n";

if ($flank == 0) {
    print STDERR "[INFO] Extract mode       : strict gene coordinates\n";
}
else {
    print STDERR "[INFO] Extract mode       : gene coordinates with +/- $flank bp flank\n";
}

# ------------------------------------------------------------
# Step 1. Read genes from merged_candidates.tsv
# ------------------------------------------------------------

my %candidate_genes;

if (-e $candidate_tsv && -s $candidate_tsv) {
    $candidate_tsv = abs_path($candidate_tsv);
    %candidate_genes = read_candidate_genes($candidate_tsv);

    my $n_candidate = scalar keys %candidate_genes;
    print STDERR "[INFO] Candidate genes found: $n_candidate\n";
}
elsif (-e $candidate_tsv) {
    print STDERR "[INFO] Candidate TSV is empty. Only HCM core genes will be extracted.\n";
}
else {
    print STDERR "[WARN] Candidate TSV not found: $candidate_tsv\n";
    print STDERR "[WARN] Only HCM core genes will be extracted.\n";
}

# ------------------------------------------------------------
# Step 2. Read HCM core genes
# ------------------------------------------------------------

my %core_genes = read_core_gene_list($core_gene_list);

my $n_core = scalar keys %core_genes;
print STDERR "[INFO] HCM core genes found: $n_core\n";

# ------------------------------------------------------------
# Step 3. Merge target genes
# ------------------------------------------------------------

my %target_genes;

for my $gene (keys %candidate_genes) {
    $target_genes{$gene} = "candidate";
}

for my $gene (keys %core_genes) {
    if (exists $target_genes{$gene}) {
        $target_genes{$gene} = "candidate,core";
    }
    else {
        $target_genes{$gene} = "core";
    }
}

my $n_target = scalar keys %target_genes;

die "[ERROR] No target gene found from candidate TSV or HCM core gene list\n"
    if $n_target == 0;

print STDERR "[INFO] Total target genes for BAM extraction: $n_target\n";

# ------------------------------------------------------------
# Step 4. Read gene coordinates and apply flank
# ------------------------------------------------------------

my %gene_regions = read_gene_regions($gene_txt, \%target_genes, $flank);

my $matched_count = scalar keys %gene_regions;

print STDERR "[INFO] Genes matched in gene TXT: $matched_count / $n_target\n";

for my $gene (sort keys %target_genes) {
    if (!exists $gene_regions{$gene}) {
        print STDERR "[WARN] Target gene has no coordinate in REFSEQ_MANE_SELECT_GENE_TXT and will be skipped: $gene\n";
    }
}

die "[ERROR] No target gene matched coordinates in $gene_txt\n"
    if $matched_count == 0;

# ------------------------------------------------------------
# Step 5. Extract BAM and generate BAM index
# ------------------------------------------------------------

my $done = 0;

for my $gene (sort keys %gene_regions) {
    my @regions = @{ $gene_regions{$gene} };
    next unless @regions;

    my $source    = $target_genes{$gene} || "unknown";
    my $safe_gene = sanitize_filename($gene);
    my $out_bam   = "$outdir/$safe_gene.bam";

    print STDERR "[INFO] Extracting gene: $gene [$source]\n";
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

    $done++;
}

print STDERR "[INFO] Gene BAM extraction finished. Total BAMs generated: $done\n";
print STDERR "[INFO] All done.\n";

exit 0;

# ============================================================
# Subroutines
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
        die "[ERROR] $name not executable: $cmd\n"
            unless -x $cmd;
    }
    else {
        my $check = system("command -v $cmd >/dev/null 2>&1");
        die "[ERROR] $name not found in PATH: $cmd\n"
            if $check != 0;
    }

    return 1;
}

sub read_candidate_genes {
    my ($file) = @_;

    open my $fh, "<", $file
        or die "[ERROR] Cannot open candidate file $file: $!\n";

    my $header = <$fh>;

    if (!defined $header) {
        close $fh;
        return ();
    }

    chomp $header;
    $header =~ s/\r$//;

    my @cols = split /\t/, $header, -1;

    my $gene_idx = find_col(
        \@cols,
        qw(
            Gene
            Annotated_Gene
            gene
            annotated_gene
            Gene_Name
            GeneName
        )
    );

    if ($gene_idx < 0) {
        print STDERR "[WARN] No Gene or Annotated_Gene column found in candidate file: $file\n";
        print STDERR "[WARN] Candidate genes will be ignored. HCM core genes will still be used.\n";

        close $fh;
        return ();
    }

    my %genes;

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;

        my @f = split /\t/, $line, -1;
        next unless defined $f[$gene_idx];

        my $gene_field = $f[$gene_idx];

        add_gene_field_to_hash($gene_field, \%genes);
    }

    close $fh;

    return %genes;
}

sub read_core_gene_list {
    my ($file) = @_;

    open my $fh, "<", $file
        or die "[ERROR] Cannot open HCM core gene list $file: $!\n";

    my %genes;
    my $line_no = 0;

    while (my $line = <$fh>) {
        $line_no++;

        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        $line =~ s/\s+#.*$//;
        $line =~ s/^\s+|\s+$//g;

        next if $line eq "";

        # HCM_CORE_GENE_LIST is expected to be a tab-delimited table:
        #   Gene    Classification
        #   PRKAG2  Definitive
        # Only the first column is used.
        my @f = split /\t/, $line, -1;
        my $gene = $f[0];

        $gene =~ s/^\s+|\s+$//g;

        next if $gene eq "";
        next if $gene eq ".";
        next if uc($gene) eq "NA";

        # Skip header line.
        next if $line_no == 1 && $gene =~ /^Gene$/i;

        $genes{$gene} = 1;
    }

    close $fh;

    return %genes;
}

sub add_gene_field_to_hash {
    my ($gene_field, $genes_ref) = @_;

    return unless defined $gene_field;

    $gene_field =~ s/^\s+|\s+$//g;

    return if $gene_field eq "";
    return if $gene_field eq ".";
    return if uc($gene_field) eq "NA";

    my @items = split /[;,]+/, $gene_field;

    for my $gene (@items) {
        $gene =~ s/^\s+|\s+$//g;

        next if $gene eq "";
        next if $gene eq ".";
        next if uc($gene) eq "NA";

        # Handle strings such as:
        #   FHOD3|NM_001281740|EX12
        #   FHOD3:NM_001281740
        # but keep normal gene symbols unchanged.
        if ($gene =~ /^([A-Za-z0-9_.-]+)\|/) {
            $gene = $1;
        }
        elsif ($gene =~ /^([A-Za-z0-9_.-]+):/) {
            $gene = $1;
        }

        $genes_ref->{$gene} = 1;
    }
}

sub read_gene_regions {
    my ($file, $target_genes_ref, $flank) = @_;

    open my $fh, "<", $file
        or die "[ERROR] Cannot open gene TXT file $file: $!\n";

    my $header = <$fh>;

    die "[ERROR] Empty gene TXT file: $file\n"
        unless defined $header;

    chomp $header;
    $header =~ s/\r$//;

    my @cols = split /\t/, $header, -1;

    my $gene_idx = find_col(
        \@cols,
        qw(Gene GeneName gene gene_name gene_symbol Symbol symbol)
    );

    my $chr_idx = find_col(
        \@cols,
        qw(Chr Chrom Chromosome chr chrom chromosome)
    );

    my $start_idx = find_col(
        \@cols,
        qw(Start Gene_start GeneStart start gene_start geneStart txStart)
    );

    my $end_idx = find_col(
        \@cols,
        qw(End Gene_end GeneEnd end gene_end geneEnd txEnd)
    );

    if ($gene_idx < 0 || $chr_idx < 0 || $start_idx < 0 || $end_idx < 0) {
        die "[ERROR] Required columns not found in gene TXT: $file\n"
          . "[ERROR] Required columns include Gene, Chrom, Start and End.\n";
    }

    my %gene_regions;

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

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

        if ($start > $end) {
            print STDERR "[WARN] Invalid region skipped for gene $gene: $chr:$start-$end\n";
            next;
        }

        my $region_start = $start - $flank;
        my $region_end   = $end + $flank;

        $region_start = 1 if $region_start < 1;

        my $region = "$chr:$region_start-$region_end";

        push @{ $gene_regions{$gene} }, $region;
    }

    close $fh;

    for my $gene (keys %gene_regions) {
        my %seen;
        my @uniq_regions;

        for my $region (@{ $gene_regions{$gene} }) {
            next if $seen{$region}++;
            push @uniq_regions, $region;
        }

        $gene_regions{$gene} = \@uniq_regions;
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
    --candidate test/test_results/25B09089386/04.candidates/25B09089386.merged_candidates.tsv \\
    --bam /path/to/sample.bam \\
    --outdir test/test_results/25B09089386/05.gene_bam \\
    --flank 0

Required:
  --config|-conf
      Config file.

  --candidate|-c
      merged_candidates.tsv file.

      Supported gene columns:
        Gene
        Annotated_Gene

      If this file is empty or has no valid gene column,
      the script will still extract genes from HCM_CORE_GENE_LIST.

  --bam|-b
      Input coordinate-sorted BAM file.

  --outdir|-o
      Output directory.

Optional:
  --flank
      Upstream/downstream flanking size in bp around gene coordinates.
      Default: 0

  --help|-h
      Show this help message.

Required config keys:
  SAMTOOLS
  REFSEQ_MANE_SELECT_GENE_TXT
  HCM_CORE_GENE_LIST

Gene extraction logic:
  target_genes = union(
      genes from merged_candidates.tsv,
      genes from HCM_CORE_GENE_LIST
  )

HCM_CORE_GENE_LIST format:
  Gene    Classification
  PRKAG2  Definitive
  ACTN2   Definitive
  FHOD3   Definitive

Only the first column is used as gene symbol.

Extraction rule:
  If --flank 0:
      strictly extract BAM by gene coordinates from REFSEQ_MANE_SELECT_GENE_TXT.

  If --flank N:
      extract gene coordinates with N bp upstream/downstream extension.

Output:
  outdir/
    GENE.bam
    GENE.bam.bai

Example:
  perl bin/extract_candidate_gene_bam.pl \\
    --config conf/hcm_exondel.example.conf \\
    --candidate test/test_results/25B09089386/04.candidates/25B09089386.merged_candidates.tsv \\
    --bam /ehpcdata/fulongfei/project/XJ_HCM_WGS_FHOD3/JX_2/25B09089386.final.merge.bam \\
    --outdir test/test_results/25B09089386/05.gene_bam \\
    --flank 0
USAGE
}

