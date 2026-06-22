#!/usr/bin/env perl
use strict;
use warnings;

# ============================================================
# Script: extract_valid_tlen.pl
#
# Purpose:
#   Extract TLEN values from valid paired-end reads in a BAM file.
#
# Usage:
#   perl extract_valid_tlen.pl <bam_abs_path> <output_abs_path>
#
# Example:
#   perl extract_valid_tlen.pl \
#       /data/WGS/sample.bam \
#       /data/WGS/sample.valid_tlen.tsv
#
# Output columns:
#   ReadName    TLEN
#
# Requirements:
#   samtools
# ============================================================

my ($bam, $out) = @ARGV;

die "Usage: perl $0 <bam_abs_path> <output_abs_path>\n"
    unless defined $bam && defined $out;

die "[ERROR] BAM file does not exist: $bam\n"
    unless -e $bam;

die "[ERROR] BAM path must be absolute: $bam\n"
    unless $bam =~ /^\//;

die "[ERROR] Output path must be absolute: $out\n"
    unless $out =~ /^\//;

# -----------------------------
# Parameters
# -----------------------------
my $MIN_MAPQ = 20;

# Require:
# 0x2 = properly paired
my $REQUIRE_FLAG = 2;

# Exclude:
# 0x4   read unmapped
# 0x8   mate unmapped
# 0x100 secondary alignment
# 0x400 duplicate
# 0x800 supplementary alignment
#
# Sum = 4 + 8 + 256 + 1024 + 2048 = 3340
my $EXCLUDE_FLAG = 3340;

# -----------------------------
# Check samtools
# -----------------------------
my $samtools = `which samtools 2>/dev/null`;
chomp $samtools;

die "[ERROR] samtools not found in PATH\n"
    unless $samtools;

# -----------------------------
# Run samtools safely
# -----------------------------
open my $SAM, "-|",
    "samtools", "view",
    "-f", $REQUIRE_FLAG,
    "-F", $EXCLUDE_FLAG,
    "-q", $MIN_MAPQ,
    $bam
    or die "[ERROR] Failed to run samtools view on BAM: $bam\n";

open my $OUT, ">", $out
    or die "[ERROR] Cannot write output file: $out\n";

# Output header
print $OUT join("\t", "ReadName", "TLEN"), "\n";

my $total_records = 0;
my $valid_tlen_count = 0;

while (my $line = <$SAM>) {
    chomp $line;
    next if $line =~ /^\s*$/;

    my @f = split /\t/, $line;

    # SAM columns, 1-based:
    # 1  QNAME
    # 2  FLAG
    # 5  MAPQ
    # 9  TLEN
    #
    # Perl array index:
    # QNAME = $f[0]
    # FLAG  = $f[1]
    # MAPQ  = $f[4]
    # TLEN  = $f[8]

    my $read_name = $f[0];
    my $tlen      = $f[8];

    $total_records++;

    next unless defined $read_name;
    next unless defined $tlen;
    next unless $tlen =~ /^-?\d+$/;

    # Use only TLEN > 0.
    # This avoids counting the same paired-end fragment twice.
    # Usually one read has positive TLEN and the mate has negative TLEN.
    next unless $tlen > 0;

    print $OUT join("\t", $read_name, $tlen), "\n";
    $valid_tlen_count++;
}

close $SAM;
close $OUT;

print "[INFO] Done.\n";
print "[INFO] BAM: $bam\n";
print "[INFO] Output: $out\n";
print "[INFO] MIN_MAPQ: $MIN_MAPQ\n";
print "[INFO] Total records after samtools filtering: $total_records\n";
print "[INFO] Valid TLEN count: $valid_tlen_count\n";

if ($valid_tlen_count == 0) {
    warn "[WARNING] No valid TLEN values were extracted. Please check BAM file and filtering criteria.\n";
}

