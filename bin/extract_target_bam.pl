#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);
use File::Basename qw(dirname);
use Cwd qw(abs_path);

# ============================================================
# Script: extract_target_bam.pl
#
# Purpose:
#   Extract one target BAM containing HCM core genes plus flank.
#
# Pipeline logic:
#   The HCMExonDel pipeline only analyzes HCM core genes.
#   This script extracts core genes + TARGET_REGION_FLANK
#   from the original WGS BAM and generates one target BAM.
#
# Design principle:
#   Command-line arguments only define runtime I/O:
#     --config
#     --sample
#     --bam
#     --outdir
#
#   All analysis parameters are read from config:
#     SAMTOOLS
#     REFSEQ_MANE_SELECT_GENE_TXT
#     HCM_CORE_GENE_LIST
#     TARGET_REGION_FLANK
#     TARGET_BAM_THREADS
#
# Coordinate convention:
#   REFSEQ_MANE_SELECT_GENE_TXT:
#     1-based closed coordinates
#
#   BED output:
#     0-based half-open coordinates
#
# Output:
#   OUTDIR/SAMPLE.target_regions.bed
#   OUTDIR/SAMPLE.target_regions.tsv
#   OUTDIR/SAMPLE.target.bam
#   OUTDIR/SAMPLE.target.bam.bai
#   OUTDIR/SAMPLE.target_bam.summary.tsv
# ============================================================

my $config;
my $sample;
my $bam;
my $outdir;
my $help = 0;

GetOptions(
    "config=s" => \$config,
    "sample=s" => \$sample,
    "bam=s"    => \$bam,
    "outdir=s" => \$outdir,
    "help"     => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage()
    unless defined $config
        && defined $sample
        && defined $bam
        && defined $outdir;

# -----------------------------
# Validate runtime I/O arguments
# -----------------------------
$config = abs_path($config) if -e $config;
die "[ERROR] Config file not found: $config\n"
    unless defined $config && -s $config;
die "[ERROR] Config path must be absolute after resolution: $config\n"
    unless $config =~ m{^/};

die "[ERROR] Invalid sample name: $sample\n"
  . "        Sample name can only contain letters, numbers, dot, underscore and hyphen.\n"
    unless defined $sample && $sample =~ /^[A-Za-z0-9_.-]+$/;

$bam = abs_path($bam) if -e $bam;
die "[ERROR] BAM file not found: $bam\n"
    unless defined $bam && -s $bam;
die "[ERROR] BAM path must be absolute after resolution: $bam\n"
    unless $bam =~ m{^/};
die "[ERROR] BAM file must end with .bam: $bam\n"
    unless $bam =~ /\.bam$/i;

die "[ERROR] Output directory path must be absolute: $outdir\n"
    unless defined $outdir && $outdir =~ m{^/};

make_path($outdir) unless -d $outdir;
$outdir = abs_path($outdir);
die "[ERROR] Cannot create output directory: $outdir\n"
    unless defined $outdir && -d $outdir;

# -----------------------------
# Read config
# -----------------------------
my %CONF = read_config($config);
my $project_root = detect_project_root($config);

my $samtools = get_conf_required(\%CONF, "SAMTOOLS");

my $gene_txt = resolve_config_path(
    get_conf_required(\%CONF, "REFSEQ_MANE_SELECT_GENE_TXT"),
    $project_root
);

my $core_gene_list = resolve_config_path(
    get_conf_required(\%CONF, "HCM_CORE_GENE_LIST"),
    $project_root
);

my $flank = get_conf_required(\%CONF, "TARGET_REGION_FLANK");
my $threads = get_conf_required(\%CONF, "TARGET_BAM_THREADS");

validate_nonnegative_integer("TARGET_REGION_FLANK", $flank);
validate_positive_integer("TARGET_BAM_THREADS", $threads);

die "[ERROR] REFSEQ_MANE_SELECT_GENE_TXT file not found: $gene_txt\n"
    unless -s $gene_txt;

die "[ERROR] HCM_CORE_GENE_LIST file not found: $core_gene_list\n"
    unless -s $core_gene_list;

check_executable($samtools, "SAMTOOLS");
check_bam_index($bam);

# -----------------------------
# Output files
# -----------------------------
my $target_bed   = "$outdir/$sample.target_regions.bed";
my $target_tsv   = "$outdir/$sample.target_regions.tsv";
my $target_bam   = "$outdir/$sample.target.bam";
my $tmp_bam      = "$outdir/$sample.target.unsorted.tmp.bam";
my $target_bai   = "$target_bam.bai";
my $summary_file = "$outdir/$sample.target_bam.summary.tsv";

# -----------------------------
# Load BAM chromosome lengths
# -----------------------------
my %bam_chr_len = read_bam_chrom_lengths($samtools, $bam);

# -----------------------------
# Load core gene list and gene annotation
# -----------------------------
my %core_gene   = read_core_gene_list($core_gene_list);
my %gene_region = read_gene_regions($gene_txt);

# -----------------------------
# Check all core genes exist in annotation
# -----------------------------
my @missing_genes;

for my $gene (sort keys %core_gene) {
    push @missing_genes, $gene unless exists $gene_region{$gene};
}

if (@missing_genes) {
    die "[ERROR] The following core genes were not found in REFSEQ_MANE_SELECT_GENE_TXT:\n"
      . join("\n", map { "        $_" } @missing_genes)
      . "\n";
}

# -----------------------------
# Generate target regions
# -----------------------------
my @raw_regions;

for my $gene (sort keys %core_gene) {
    my $r = $gene_region{$gene};

    my $chr = $r->{chrom};

    die "[ERROR] Chromosome '$chr' for gene '$gene' is not found in BAM header/index.\n"
      . "        Please check whether annotation chromosome names match the BAM.\n"
        unless exists $bam_chr_len{$chr};

    my $chr_len = $bam_chr_len{$chr};

    my $target_start = $r->{start} - $flank;
    my $target_end   = $r->{end}   + $flank;

    $target_start = 1 if $target_start < 1;
    $target_end = $chr_len if $target_end > $chr_len;

    die "[ERROR] Invalid target interval after flank clipping for gene $gene: $chr:$target_start-$target_end\n"
        unless $target_start <= $target_end;

    push @raw_regions, {
        gene         => $gene,
        transcript   => $r->{transcript},
        chrom        => $chr,
        gene_start   => $r->{start},
        gene_end     => $r->{end},
        target_start => $target_start,
        target_end   => $target_end,
        flank        => $flank,
    };
}

die "[ERROR] No target regions generated. Please check HCM_CORE_GENE_LIST and REFSEQ_MANE_SELECT_GENE_TXT.\n"
    unless @raw_regions;

my @merged_regions = merge_regions(\@raw_regions);

write_region_tsv($target_tsv, \@raw_regions);
write_bed($target_bed, \@merged_regions);

# -----------------------------
# Extract target BAM
# -----------------------------
unlink $tmp_bam if -e $tmp_bam;
unlink $target_bam if -e $target_bam;
unlink $target_bai if -e $target_bai;
unlink "$target_bam.csi" if -e "$target_bam.csi";

run_cmd(
    $samtools, "view",
    "-@", $threads,
    "-bh",
    "-M",
    "-L", $target_bed,
    "-o", $tmp_bam,
    $bam
);

die "[ERROR] Temporary target BAM was not generated: $tmp_bam\n"
    unless -s $tmp_bam;

run_cmd(
    $samtools, "sort",
    "-@", $threads,
    "-o", $target_bam,
    $tmp_bam
);

die "[ERROR] Target BAM was not generated: $target_bam\n"
    unless -s $target_bam;

run_cmd(
    $samtools, "index",
    "-@", $threads,
    $target_bam
);

die "[ERROR] Target BAM index was not generated: $target_bai\n"
    unless -s $target_bai;

unlink $tmp_bam if -e $tmp_bam;

write_summary(
    $summary_file,
    sample                => $sample,
    config                => $config,
    input_bam             => $bam,
    output_bam            => $target_bam,
    output_bam_index      => $target_bai,
    target_bed            => $target_bed,
    target_region_tsv     => $target_tsv,
    core_gene_list        => $core_gene_list,
    gene_annotation       => $gene_txt,
    target_region_flank   => $flank,
    target_bam_threads    => $threads,
    n_core_genes          => scalar(keys %core_gene),
    n_raw_gene_regions    => scalar(@raw_regions),
    n_merged_bed_regions  => scalar(@merged_regions),
);

print "[INFO] Target BAM extraction finished.\n";
print "[INFO] Sample: $sample\n";
print "[INFO] Input BAM: $bam\n";
print "[INFO] Output BAM: $target_bam\n";
print "[INFO] Output BAM index: $target_bai\n";
print "[INFO] Target BED: $target_bed\n";
print "[INFO] Target region TSV: $target_tsv\n";
print "[INFO] Summary: $summary_file\n";
print "[INFO] TARGET_REGION_FLANK: $flank\n";
print "[INFO] TARGET_BAM_THREADS: $threads\n";

exit 0;

# ============================================================
# Config helpers
# ============================================================

sub read_config {
    my ($file) = @_;

    my %conf;

    open my $FH, "<", $file
        or die "[ERROR] Cannot open config file: $file\n";

    while (my $line = <$FH>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        # Remove inline comments:
        # KEY=value   # comment
        $line =~ s/\s+#.*$//;

        next unless $line =~ /^\s*([^=\s]+)\s*=\s*(.*?)\s*$/;

        my ($key, $val) = ($1, $2);

        $key =~ s/^\s+|\s+$//g;
        $val =~ s/^\s+|\s+$//g;

        # Remove one pair of simple quotes if present
        $val =~ s/^['"]//;
        $val =~ s/['"]$//;

        die "[ERROR] Empty config key found in $file\n"
            if $key eq "";

        $conf{$key} = $val;
    }

    close $FH;

    return %conf;
}

sub get_conf_required {
    my ($conf_ref, $key) = @_;

    die "[ERROR] Required config parameter missing: $key\n"
        unless exists $conf_ref->{$key}
            && defined $conf_ref->{$key}
            && $conf_ref->{$key} ne "";

    return $conf_ref->{$key};
}

sub detect_project_root {
    my ($config_path) = @_;

    my $config_dir = dirname($config_path);

    if ($config_dir =~ m{/conf$}) {
        return abs_path(dirname($config_dir));
    }

    return abs_path($config_dir);
}

sub resolve_config_path {
    my ($path, $project_root) = @_;

    die "[ERROR] Empty config path provided\n"
        unless defined $path && $path ne "";

    if ($path =~ m{^/}) {
        my $abs = abs_path($path);

        die "[ERROR] Path not found: $path\n"
            unless defined $abs;

        return $abs;
    }

    my $full = "$project_root/$path";
    my $abs = abs_path($full);

    die "[ERROR] Path not found: $full\n"
        unless defined $abs;

    return $abs;
}

# ============================================================
# Validation helpers
# ============================================================

sub validate_positive_integer {
    my ($name, $v) = @_;

    die "[ERROR] $name must be a positive integer. Observed: $v\n"
        unless defined $v && $v =~ /^\d+$/ && $v >= 1;

    return 1;
}

sub validate_nonnegative_integer {
    my ($name, $v) = @_;

    die "[ERROR] $name must be a non-negative integer. Observed: $v\n"
        unless defined $v && $v =~ /^\d+$/ && $v >= 0;

    return 1;
}

sub check_executable {
    my ($cmd, $name) = @_;

    die "[ERROR] Empty executable path for $name\n"
        unless defined $cmd && $cmd ne "";

    if ($cmd =~ m{/}) {
        die "[ERROR] $name executable not found or not executable: $cmd\n"
            unless -x $cmd;
    } else {
        my $ret = system("command -v $cmd >/dev/null 2>&1");

        die "[ERROR] $name executable not found in PATH: $cmd\n"
            if $ret != 0;
    }

    return 1;
}

sub check_bam_index {
    my ($bam_file) = @_;

    my $bai1 = "$bam_file.bai";

    my $bai2 = $bam_file;
    $bai2 =~ s/\.bam$/.bai/i;

    die "[ERROR] BAM index not found for: $bam_file\n"
      . "        Expected one of:\n"
      . "        $bai1\n"
      . "        $bai2\n"
        unless -s $bai1 || -s $bai2;

    return 1;
}

# ============================================================
# Input parsers
# ============================================================

sub read_core_gene_list {
    my ($file) = @_;

    my %gene;

    open my $FH, "<", $file
        or die "[ERROR] Cannot open HCM_CORE_GENE_LIST: $file\n";

    my $line_no = 0;

    while (my $line = <$FH>) {
        chomp $line;
        $line =~ s/\r$//;
        $line_no++;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        my @f = split /\t/, $line;
        my $g = $f[0];

        $g = "" unless defined $g;
        $g =~ s/^\s+|\s+$//g;

        die "[ERROR] Empty gene symbol in HCM_CORE_GENE_LIST at line $line_no\n"
            if $g eq "";

        die "[ERROR] Invalid gene symbol in HCM_CORE_GENE_LIST at line $line_no: $g\n"
            unless $g =~ /^[A-Za-z0-9_.-]+$/;

        $gene{$g} = 1;
    }

    close $FH;

    die "[ERROR] No core genes found in HCM_CORE_GENE_LIST: $file\n"
        unless keys %gene;

    return %gene;
}

sub read_gene_regions {
    my ($file) = @_;

    open my $FH, "<", $file
        or die "[ERROR] Cannot open REFSEQ_MANE_SELECT_GENE_TXT: $file\n";

    my $header = <$FH>;

    die "[ERROR] Empty REFSEQ_MANE_SELECT_GENE_TXT: $file\n"
        unless defined $header;

    chomp $header;
    $header =~ s/\r$//;

    my @cols = split /\t/, $header, -1;
    my %idx;

    for my $i (0 .. $#cols) {
        my $col = $cols[$i];
        $col =~ s/^\s+|\s+$//g;
        $idx{$col} = $i;
    }

    for my $required (qw/Gene Transcript Chrom Start End Strand/) {
        die "[ERROR] Required column '$required' not found in REFSEQ_MANE_SELECT_GENE_TXT: $file\n"
            unless exists $idx{$required};
    }

    my %region;
    my $line_no = 1;

    while (my $line = <$FH>) {
        chomp $line;
        $line =~ s/\r$//;
        $line_no++;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        my @f = split /\t/, $line, -1;

        my $gene = $f[$idx{Gene}];
        my $tx   = $f[$idx{Transcript}];
        my $chr  = $f[$idx{Chrom}];
        my $s    = $f[$idx{Start}];
        my $e    = $f[$idx{End}];

        for ($gene, $tx, $chr, $s, $e) {
            $_ = "" unless defined $_;
            s/^\s+|\s+$//g;
        }

        die "[ERROR] Invalid gene annotation at line $line_no: empty Gene/Transcript/Chrom/Start/End\n"
            unless $gene ne ""
                && $tx ne ""
                && $chr ne ""
                && $s ne ""
                && $e ne "";

        die "[ERROR] Invalid Start at line $line_no: $s\n"
            unless $s =~ /^\d+$/ && $s >= 1;

        die "[ERROR] Invalid End at line $line_no: $e\n"
            unless $e =~ /^\d+$/ && $e >= 1;

        die "[ERROR] Start > End at line $line_no: $chr:$s-$e\n"
            if $s > $e;

        if (!exists $region{$gene}) {
            $region{$gene} = {
                chrom      => $chr,
                start      => $s,
                end        => $e,
                transcript => $tx,
            };
        } else {
            die "[ERROR] Gene $gene appears on multiple chromosomes in REFSEQ_MANE_SELECT_GENE_TXT.\n"
                if $region{$gene}->{chrom} ne $chr;

            $region{$gene}->{start} = $s
                if $s < $region{$gene}->{start};

            $region{$gene}->{end} = $e
                if $e > $region{$gene}->{end};

            $region{$gene}->{transcript} .= ",$tx"
                unless $region{$gene}->{transcript} =~ /(?:^|,)\Q$tx\E(?:,|$)/;
        }
    }

    close $FH;

    die "[ERROR] No gene records found in REFSEQ_MANE_SELECT_GENE_TXT: $file\n"
        unless keys %region;

    return %region;
}

sub read_bam_chrom_lengths {
    my ($samtools, $bam_file) = @_;

    my %len;

    open my $IDX, "-|", $samtools, "idxstats", $bam_file
        or die "[ERROR] Failed to run samtools idxstats on BAM: $bam_file\n";

    while (my $line = <$IDX>) {
        chomp $line;

        my @f = split /\t/, $line;
        next unless @f >= 2;

        my ($chr, $length) = @f[0, 1];

        next if !defined $chr || $chr eq "*";
        next unless defined $length && $length =~ /^\d+$/ && $length > 0;

        $len{$chr} = $length;
    }

    close $IDX
        or die "[ERROR] samtools idxstats failed for BAM: $bam_file\n";

    die "[ERROR] No chromosome length information was obtained from BAM: $bam_file\n"
        unless keys %len;

    return %len;
}

# ============================================================
# Region processing
# ============================================================

sub chrom_order_key {
    my ($chr) = @_;

    my $x = $chr;
    $x =~ s/^chr//i;

    return sprintf("%03d", $x) if $x =~ /^\d+$/;

    my $u = uc($x);

    return "023" if $u eq "X";
    return "024" if $u eq "Y";
    return "025" if $u eq "M" || $u eq "MT";

    return "999_$chr";
}

sub merge_regions {
    my ($regions_ref) = @_;

    my @sorted = sort {
        chrom_order_key($a->{chrom}) cmp chrom_order_key($b->{chrom})
            || $a->{target_start} <=> $b->{target_start}
            || $a->{target_end} <=> $b->{target_end}
            || $a->{gene} cmp $b->{gene}
    } @$regions_ref;

    my @merged;

    for my $r (@sorted) {
        if (
            !@merged
            || $merged[-1]->{chrom} ne $r->{chrom}
            || $r->{target_start} > $merged[-1]->{end} + 1
        ) {
            push @merged, {
                chrom => $r->{chrom},
                start => $r->{target_start},
                end   => $r->{target_end},
                genes => [ $r->{gene} ],
            };
        } else {
            $merged[-1]->{end} = $r->{target_end}
                if $r->{target_end} > $merged[-1]->{end};

            push @{ $merged[-1]->{genes} }, $r->{gene};
        }
    }

    return @merged;
}

sub write_region_tsv {
    my ($file, $regions_ref) = @_;

    open my $OUT, ">", $file
        or die "[ERROR] Cannot write target region TSV: $file\n";

    print $OUT join(
        "\t",
        qw/
          Gene
          Transcript
          Chrom
          Gene_Start
          Gene_End
          Flank
          Target_Start
          Target_End
          Target_Length
        /
    ), "\n";

    for my $r (
        sort {
            chrom_order_key($a->{chrom}) cmp chrom_order_key($b->{chrom})
                || $a->{target_start} <=> $b->{target_start}
                || $a->{target_end} <=> $b->{target_end}
                || $a->{gene} cmp $b->{gene}
        } @$regions_ref
    ) {
        my $len = $r->{target_end} - $r->{target_start} + 1;

        print $OUT join(
            "\t",
            $r->{gene},
            $r->{transcript},
            $r->{chrom},
            $r->{gene_start},
            $r->{gene_end},
            $r->{flank},
            $r->{target_start},
            $r->{target_end},
            $len,
        ), "\n";
    }

    close $OUT;
}

sub write_bed {
    my ($file, $regions_ref) = @_;

    open my $BED, ">", $file
        or die "[ERROR] Cannot write target BED: $file\n";

    for my $r (@$regions_ref) {
        my $bed_start = $r->{start} - 1;
        my $bed_end   = $r->{end};
        my $name      = join(",", @{ $r->{genes} });

        die "[ERROR] Invalid BED interval generated: $r->{chrom}:$bed_start-$bed_end\n"
            unless $bed_start >= 0 && $bed_start < $bed_end;

        print $BED join(
            "\t",
            $r->{chrom},
            $bed_start,
            $bed_end,
            $name
        ), "\n";
    }

    close $BED;
}

# ============================================================
# Command and summary helpers
# ============================================================

sub run_cmd {
    my @cmd = @_;

    print "[CMD] ", join(" ", map { shell_quote($_) } @cmd), "\n";

    system(@cmd) == 0
        or die "[ERROR] Command failed: "
        . join(" ", map { shell_quote($_) } @cmd)
        . "\n";
}

sub shell_quote {
    my ($str) = @_;

    die "[ERROR] Undefined shell argument\n"
        unless defined $str;

    return "''" if $str eq "";

    return $str
        if $str =~ /^[A-Za-z0-9_\.\-\/\:=,\+]+$/;

    $str =~ s/'/'\\''/g;

    return "'$str'";
}

sub write_summary {
    my ($file, %kv) = @_;

    open my $OUT, ">", $file
        or die "[ERROR] Cannot write summary file: $file\n";

    print $OUT "Key\tValue\n";

    for my $k (sort keys %kv) {
        print $OUT "$k\t$kv{$k}\n";
    }

    close $OUT;
}

sub usage {
    return <<'USAGE';
Usage:
  perl bin/extract_target_bam.pl \
    --config /abs/path/conf/hcm_exondel.conf \
    --sample SAMPLE_ID \
    --bam /abs/path/sample.wgs.bam \
    --outdir /abs/path/results/SAMPLE/00.target_bam

Required config parameters:
  SAMTOOLS
  REFSEQ_MANE_SELECT_GENE_TXT
  HCM_CORE_GENE_LIST
  TARGET_REGION_FLANK
  TARGET_BAM_THREADS

Notes:
  1. The HCMExonDel pipeline only analyzes HCM core genes.
  2. This script extracts HCM core genes plus TARGET_REGION_FLANK into one target BAM.
  3. All analysis parameters are read from the config file.
  4. Command-line arguments only define runtime input/output identity.
  5. The output target BAM is sorted and indexed.
USAGE
}

