#!/usr/bin/env perl
# -*- coding: utf-8 -*-

use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(dirname);
use File::Path qw(make_path);

# ============================================================
# check_bam.pl
#
# Check BAM file before running HCMExonDel.
#
# Required:
#   --config
#   --bam
#
# Optional:
#   --sample
#   --out
#
# Config parameters:
#   SAMTOOLS
#   TARGET_REGION_BED
#   CHECK_TARGET_READS
# ============================================================

my $config;
my $bam;
my $sample = "NA";
my $out;
my $help = 0;

GetOptions(
    "config=s" => \$config,
    "bam=s"    => \$bam,
    "sample=s" => \$sample,
    "out=s"    => \$out,
    "help"     => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless $config && $bam;

$config = abs_path($config);
$bam    = abs_path($bam);

die "[ERROR] Config file not found: $config\n" unless -s $config;
die "[ERROR] BAM file not found: $bam\n" unless -s $bam;

my %CONF = read_config($config);

my $samtools = get_conf_required(\%CONF, "SAMTOOLS");
my $target_bed = get_conf_value(\%CONF, "TARGET_REGION_BED", "");
my $check_target_reads = normalize_bool(
    get_conf_value(\%CONF, "CHECK_TARGET_READS", 1)
);

if ($target_bed) {
    $target_bed = abs_path($target_bed);
    die "[ERROR] TARGET_REGION_BED file not found: $target_bed\n" unless -s $target_bed;
}

if (!$out) {
    $out = $sample ne "NA" ? "$sample.bam_check.tsv" : "bam_check.tsv";
}

my $outdir = dirname($out);
make_path($outdir) unless -d $outdir;

my @results;

# ============================================================
# Check 1. samtools
# ============================================================

my $samtools_status = check_samtools($samtools);
push @results, result_row(
    "samtools_available",
    $samtools_status ? "PASS" : "FAIL",
    $samtools
);

if (!$samtools_status) {
    write_report($out, \@results);
    die "[ERROR] samtools is not available: $samtools\n";
}

# ============================================================
# Check 2. BAM existence
# ============================================================

push @results, result_row(
    "bam_exists",
    -s $bam ? "PASS" : "FAIL",
    $bam
);

# ============================================================
# Check 3. BAM index
# ============================================================

my $bai = find_bam_index($bam);

push @results, result_row(
    "bam_index_exists",
    $bai ? "PASS" : "FAIL",
    $bai || "BAM index not found"
);

# ============================================================
# Check 4. BAM header readable
# ============================================================

my $header = get_bam_header($samtools, $bam);
my $header_ok = $header ? 1 : 0;

push @results, result_row(
    "bam_header_readable",
    $header_ok ? "PASS" : "FAIL",
    $header_ok ? "BAM header can be read" : "Failed to read BAM header"
);

if (!$header_ok) {
    write_report($out, \@results);
    die "[ERROR] Failed to read BAM header: $bam\n";
}

# ============================================================
# Check 5. Coordinate sorted
# ============================================================

my $sort_order = parse_sort_order($header);

push @results, result_row(
    "bam_sort_order",
    $sort_order eq "coordinate" ? "PASS" : "FAIL",
    "SO=$sort_order"
);

# ============================================================
# Check 6. Chromosome names in BAM
# ============================================================

my %bam_chrom = parse_bam_chromosomes($header);

push @results, result_row(
    "bam_chromosome_count",
    scalar(keys %bam_chrom) > 0 ? "PASS" : "FAIL",
    scalar(keys %bam_chrom) . " chromosomes found in BAM header"
);

# ============================================================
# Check 7. TARGET_REGION_BED chromosome consistency
# ============================================================

if ($target_bed) {
    my %bed_chrom = parse_bed_chromosomes($target_bed);

    my @missing;
    foreach my $chr (sort keys %bed_chrom) {
        push @missing, $chr unless exists $bam_chrom{$chr};
    }

    if (@missing) {
        push @results, result_row(
            "target_bed_chromosome_consistency",
            "FAIL",
            "Chromosomes in TARGET_REGION_BED not found in BAM: " . join(",", @missing)
        );
    }
    else {
        push @results, result_row(
            "target_bed_chromosome_consistency",
            "PASS",
            "All TARGET_REGION_BED chromosomes are found in BAM header"
        );
    }
}
else {
    push @results, result_row(
        "target_bed_chromosome_consistency",
        "SKIP",
        "TARGET_REGION_BED not provided in config"
    );
}

# ============================================================
# Check 8. Optional target region read check
# ============================================================

if ($target_bed && $check_target_reads) {
    my ($regions_checked, $regions_with_reads) =
        check_target_region_reads($samtools, $bam, $target_bed);

    my $status = $regions_with_reads > 0 ? "PASS" : "WARN";

    push @results, result_row(
        "target_region_reads",
        $status,
        "Regions checked=$regions_checked; regions with reads=$regions_with_reads"
    );
}
else {
    push @results, result_row(
        "target_region_reads",
        "SKIP",
        "CHECK_TARGET_READS=0 or TARGET_REGION_BED not provided"
    );
}

# ============================================================
# Final status
# ============================================================

my $final_status = "PASS";

foreach my $r (@results) {
    if ($r->{Status} eq "FAIL") {
        $final_status = "FAIL";
        last;
    }
}

push @results, result_row(
    "final_status",
    $final_status,
    "BAM check completed"
);

write_report($out, \@results);

print "[INFO] BAM check finished\n";
print "[INFO] Sample       : $sample\n";
print "[INFO] BAM          : $bam\n";
print "[INFO] Report       : $out\n";
print "[INFO] Final status : $final_status\n";

if ($final_status eq "FAIL") {
    exit 1;
}

exit 0;


# ============================================================
# Subroutines
# ============================================================

sub check_samtools {
    my ($samtools) = @_;

    my $cmd = shell_quote($samtools) . " --version 2>/dev/null";
    my $ret = system($cmd);

    return $ret == 0 ? 1 : 0;
}


sub find_bam_index {
    my ($bam) = @_;

    my $bai1 = "$bam.bai";

    my $bai2 = $bam;
    $bai2 =~ s/\.bam$/.bai/;

    my $csi = "$bam.csi";

    return $bai1 if -s $bai1;
    return $bai2 if -s $bai2;
    return $csi  if -s $csi;

    return "";
}


sub get_bam_header {
    my ($samtools, $bam) = @_;

    my $cmd = join(" ",
        shell_quote($samtools),
        "view",
        "-H",
        shell_quote($bam)
    );

    my $header = `$cmd 2>/dev/null`;

    return $header;
}


sub parse_sort_order {
    my ($header) = @_;

    foreach my $line (split /\n/, $header) {
        next unless $line =~ /^\@HD/;

        if ($line =~ /\bSO:([^\t]+)/) {
            return $1;
        }
    }

    return "unknown";
}


sub parse_bam_chromosomes {
    my ($header) = @_;

    my %chrom;

    foreach my $line (split /\n/, $header) {
        next unless $line =~ /^\@SQ/;

        my $chr;
        my $len;

        if ($line =~ /\bSN:([^\t]+)/) {
            $chr = $1;
        }

        if ($line =~ /\bLN:(\d+)/) {
            $len = $1;
        }

        if ($chr) {
            $chrom{$chr} = $len || 0;
        }
    }

    return %chrom;
}


sub parse_bed_chromosomes {
    my ($bed) = @_;

    my %chrom;

    open my $fh, "<", $bed
        or die "[ERROR] Cannot open TARGET_REGION_BED: $bed\n";

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        my @f = split /\t/, $line;

        next unless @f >= 3;

        my $chr = $f[0];

        $chrom{$chr} = 1 if $chr;
    }

    close $fh;

    return %chrom;
}


sub check_target_region_reads {
    my ($samtools, $bam, $bed) = @_;

    my $regions_checked = 0;
    my $regions_with_reads = 0;

    open my $fh, "<", $bed
        or die "[ERROR] Cannot open TARGET_REGION_BED: $bed\n";

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        my @f = split /\t/, $line;
        next unless @f >= 3;

        my ($chr, $start0, $end1) = @f[0, 1, 2];

        next unless $start0 =~ /^\d+$/ && $end1 =~ /^\d+$/ && $end1 > $start0;

        my $start1 = $start0 + 1;
        my $region = "$chr:$start1-$end1";

        $regions_checked++;

        my $cmd = join(" ",
            shell_quote($samtools),
            "view",
            "-c",
            shell_quote($bam),
            shell_quote($region)
        );

        my $count = `$cmd 2>/dev/null`;
        chomp $count;

        if (defined $count && $count =~ /^\d+$/ && $count > 0) {
            $regions_with_reads++;
        }
    }

    close $fh;

    return ($regions_checked, $regions_with_reads);
}


sub result_row {
    my ($item, $status, $message) = @_;

    return {
        Item    => $item,
        Status  => $status,
        Message => $message,
    };
}


sub write_report {
    my ($out, $results_ref) = @_;

    open my $fh, ">", $out
        or die "[ERROR] Cannot write report: $out\n";

    print $fh join("\t", qw(SampleID Item Status Message)) . "\n";

    foreach my $r (@$results_ref) {
        print $fh join("\t",
            $sample,
            $r->{Item},
            $r->{Status},
            $r->{Message},
        ) . "\n";
    }

    close $fh;
}


sub read_config {
    my ($file) = @_;

    my %conf;

    open my $fh, "<", $file
        or die "[ERROR] Cannot open config file: $file\n";

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


sub get_conf_required {
    my ($conf_ref, $key) = @_;

    die "[ERROR] Required config parameter missing: $key\n"
        unless exists $conf_ref->{$key} && $conf_ref->{$key} ne "";

    return $conf_ref->{$key};
}


sub get_conf_value {
    my ($conf_ref, $key, $default) = @_;

    if (exists $conf_ref->{$key} && $conf_ref->{$key} ne "") {
        return $conf_ref->{$key};
    }

    return $default;
}


sub normalize_bool {
    my ($value) = @_;

    return 0 unless defined $value;

    return 1 if $value =~ /^1$/;
    return 1 if $value =~ /^true$/i;
    return 1 if $value =~ /^yes$/i;
    return 1 if $value =~ /^y$/i;

    return 0;
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

  perl scripts/check_bam.pl \\
    --config conf/hcm_exondel.conf \\
    --bam sample.sorted.bam \\
    --sample SAMPLE001 \\
    --out results/SAMPLE001/00.log/SAMPLE001.bam_check.tsv

Required arguments:

  --config       Config file
  --bam          BAM file to check

Optional arguments:

  --sample       Sample ID [default: NA]
  --out          Output check report [default: SAMPLE.bam_check.tsv]
  --help         Show this help message

Config parameters used:

  SAMTOOLS
  TARGET_REGION_BED
  CHECK_TARGET_READS

Checks performed:

  1. samtools availability
  2. BAM existence
  3. BAM index existence
  4. BAM header readability
  5. BAM sort order
  6. chromosome consistency between BAM and TARGET_REGION_BED
  7. optional target-region read count check

Recommended config:

  SAMTOOLS=/path/to/samtools
  TARGET_REGION_BED=/path/to/db/target_region.bed
  CHECK_TARGET_READS=1

USAGE
}

