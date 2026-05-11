#!/usr/bin/env perl
# -*- coding: utf-8 -*-

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path);

# ============================================================
# run_discordant_reads.pl
#
# Discordant read-pair screening from coordinate-sorted WGS BAM.
#
# Design principle:
#   Command line only provides sample-level input/output.
#   All software paths, database paths and thresholds are read
#   from config file.
#
# Required command-line arguments:
#   --config
#   --bam
#   --sample
#   --out
#
# Config-controlled parameters:
#   SAMTOOLS
#   TARGET_REGION_BED
#   SCAN_WHOLE_BAM
#   MIN_MAPQ
#   MIN_DISCORDANT_INSERT_SIZE
#   MIN_DISCORDANT_READS
#   DISCORDANT_CLUSTER_DISTANCE
#   MIN_DISCORDANT_CANDIDATE_SIZE
#   MAX_DISCORDANT_CANDIDATE_SIZE
#   EXCLUDE_DUPLICATES
# ============================================================

my $config;
my $bam;
my $sample;
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

die usage() unless $config && $bam && $sample && $out;

$config = abs_path($config);
$bam    = abs_path($bam);

die "[ERROR] Config file not found: $config\n" unless -s $config;
die "[ERROR] BAM file not found: $bam\n" unless -s $bam;

my %CONF = read_config($config);

# ============================================================
# Read all analysis parameters from config
# ============================================================

my $samtools = get_conf_required(\%CONF, "SAMTOOLS");

my $scan_whole_bam = get_conf_value(\%CONF, "SCAN_WHOLE_BAM", 0);
$scan_whole_bam = normalize_bool($scan_whole_bam);

my $target_bed = get_conf_value(\%CONF, "TARGET_REGION_BED", "");

if (!$scan_whole_bam) {
    die "[ERROR] TARGET_REGION_BED is required in config when SCAN_WHOLE_BAM=0\n"
        unless $target_bed;

    $target_bed = abs_path($target_bed);
    die "[ERROR] TARGET_REGION_BED file not found: $target_bed\n"
        unless -s $target_bed;
}

my $min_mapq = get_conf_value(\%CONF, "MIN_MAPQ", 20);

my $min_discordant_insert_size =
    get_conf_value(\%CONF, "MIN_DISCORDANT_INSERT_SIZE", 1000);

my $min_discordant_reads =
    get_conf_value(\%CONF, "MIN_DISCORDANT_READS", 3);

my $cluster_distance =
    get_conf_value(\%CONF, "DISCORDANT_CLUSTER_DISTANCE", 1000);

my $min_candidate_size =
    get_conf_value(\%CONF, "MIN_DISCORDANT_CANDIDATE_SIZE", 50);

my $max_candidate_size =
    get_conf_value(\%CONF, "MAX_DISCORDANT_CANDIDATE_SIZE", 1000000);

my $exclude_duplicates =
    normalize_bool(get_conf_value(\%CONF, "EXCLUDE_DUPLICATES", 1));

# ============================================================
# Prepare output
# ============================================================

my $outdir = dirname($out);
make_path($outdir) unless -d $outdir;

my $prefix = $out;
$prefix =~ s/\.tsv$//;

my $supporting_pairs_out = $prefix . ".supporting_pairs.tsv";
my $raw_pairs_out        = $prefix . ".raw_discordant_pairs.tsv";

# ============================================================
# Prepare scan regions
# ============================================================

my @scan_regions;

if ($scan_whole_bam) {
    push @scan_regions, {
        chr    => "ALL",
        start  => 0,
        end    => 0,
        name   => "whole_bam",
        region => "",
    };
}
else {
    @scan_regions = read_target_bed($target_bed);

    die "[ERROR] No valid target regions found in TARGET_REGION_BED: $target_bed\n"
        unless @scan_regions;
}

print "[INFO] Discordant read screening started\n";
print "[INFO] Sample                    : $sample\n";
print "[INFO] BAM                       : $bam\n";
print "[INFO] Config                    : $config\n";
print "[INFO] Scan whole BAM             : $scan_whole_bam\n";
print "[INFO] Target BED                 : " . ($target_bed || "NA") . "\n";
print "[INFO] Scan region number         : " . scalar(@scan_regions) . "\n";
print "[INFO] Min MAPQ                   : $min_mapq\n";
print "[INFO] Min discordant insert size : $min_discordant_insert_size\n";
print "[INFO] Min discordant reads       : $min_discordant_reads\n";
print "[INFO] Cluster distance           : $cluster_distance\n";

# ============================================================
# Scan BAM for discordant pairs
# ============================================================

my @discordant_pairs;

open my $raw_fh, ">", $raw_pairs_out
    or die "[ERROR] Cannot write raw pairs output: $raw_pairs_out\n";

print $raw_fh join("\t",
    qw(
        SampleID
        Read_Name
        Chrom
        Read_Pos
        Mate_Chrom
        Mate_Pos
        Insert_Size
        MAPQ
        FLAG
        CIGAR
        Target_Name
        Pair_Key
    )
) . "\n";

foreach my $region (@scan_regions) {
    my @pairs = scan_region_for_discordant_pairs(
        bam             => $bam,
        samtools        => $samtools,
        sample          => $sample,
        region          => $region,
        min_mapq        => $min_mapq,
        min_insert_size => $min_discordant_insert_size,
        exclude_dup     => $exclude_duplicates,
        raw_fh          => $raw_fh,
    );

    push @discordant_pairs, @pairs;
}

close $raw_fh;

my @unique_pairs = remove_duplicate_pairs(@discordant_pairs);

# ============================================================
# Cluster discordant pairs
# ============================================================

my @clusters = cluster_discordant_pairs(
    pairs            => \@unique_pairs,
    cluster_distance => $cluster_distance,
);

# ============================================================
# Output summary and supporting reads
# ============================================================

open my $out_fh, ">", $out
    or die "[ERROR] Cannot write output file: $out\n";

open my $support_fh, ">", $supporting_pairs_out
    or die "[ERROR] Cannot write supporting pairs output: $supporting_pairs_out\n";

print $out_fh join("\t",
    qw(
        SampleID
        Chrom
        Start
        End
        Candidate_Length
        Discordant_Reads
        Median_Insert_Size
        Min_Insert_Size
        Max_Insert_Size
        Left_Boundary_Min
        Left_Boundary_Max
        Right_Boundary_Min
        Right_Boundary_Max
        Target_Name
        Discordant_Status
        Candidate_Status
        Comment
    )
) . "\n";

print $support_fh join("\t",
    qw(
        SampleID
        Cluster_ID
        Read_Name
        Chrom
        Read_Pos
        Mate_Chrom
        Mate_Pos
        Insert_Size
        MAPQ
        FLAG
        CIGAR
        Target_Name
    )
) . "\n";

my $cluster_id = 0;
my $pass_cluster_count = 0;

foreach my $cluster (@clusters) {
    next if $cluster->{discordant_reads} < $min_discordant_reads;
    next if $cluster->{candidate_length} < $min_candidate_size;
    next if $cluster->{candidate_length} > $max_candidate_size;

    $cluster_id++;
    $pass_cluster_count++;

    my $status = "Discordant_supported";
    my $candidate_status = "Discordant_candidate";
    my $comment = "Candidate deletion region supported by discordant read pairs";

    print $out_fh join("\t",
        $sample,
        $cluster->{chr},
        $cluster->{start},
        $cluster->{end},
        $cluster->{candidate_length},
        $cluster->{discordant_reads},
        sprintf("%.2f", $cluster->{median_insert_size}),
        $cluster->{min_insert_size},
        $cluster->{max_insert_size},
        $cluster->{left_min},
        $cluster->{left_max},
        $cluster->{right_min},
        $cluster->{right_max},
        $cluster->{target_name},
        $status,
        $candidate_status,
        $comment,
    ) . "\n";

    foreach my $p (@{ $cluster->{pairs} }) {
        print $support_fh join("\t",
            $sample,
            "DR_CLUSTER_" . $cluster_id,
            $p->{read_name},
            $p->{chr},
            $p->{read_pos},
            $p->{mate_chr},
            $p->{mate_pos},
            $p->{insert_size},
            $p->{mapq},
            $p->{flag},
            $p->{cigar},
            $p->{target_name},
        ) . "\n";
    }
}

close $out_fh;
close $support_fh;

print "[INFO] Discordant reads screening finished\n";
print "[INFO] Raw discordant pairs output : $raw_pairs_out\n";
print "[INFO] Summary output              : $out\n";
print "[INFO] Supporting pairs output     : $supporting_pairs_out\n";
print "[INFO] Raw discordant pairs        : " . scalar(@discordant_pairs) . "\n";
print "[INFO] Unique discordant pairs     : " . scalar(@unique_pairs) . "\n";
print "[INFO] Total clusters              : " . scalar(@clusters) . "\n";
print "[INFO] Passed clusters             : $pass_cluster_count\n";

exit 0;


# ============================================================
# Subroutines
# ============================================================

sub scan_region_for_discordant_pairs {
    my %args = @_;

    my $bam             = $args{bam};
    my $samtools        = $args{samtools};
    my $sample          = $args{sample};
    my $region          = $args{region};
    my $min_mapq        = $args{min_mapq};
    my $min_insert_size = $args{min_insert_size};
    my $exclude_dup     = $args{exclude_dup};
    my $raw_fh          = $args{raw_fh};

    my @pairs;

    my $cmd;

    if ($region->{region}) {
        $cmd = join(" ",
            shell_quote($samtools),
            "view",
            shell_quote($bam),
            shell_quote($region->{region})
        );
    }
    else {
        $cmd = join(" ",
            shell_quote($samtools),
            "view",
            shell_quote($bam)
        );
    }

    open my $pipe, "$cmd |"
        or die "[ERROR] Failed to run command: $cmd\n";

    while (my $line = <$pipe>) {
        chomp $line;
        next if $line =~ /^\s*$/;

        my @f = split /\t/, $line;
        next unless @f >= 11;

        my $qname = $f[0];
        my $flag  = $f[1];
        my $rname = $f[2];
        my $pos   = $f[3];
        my $mapq  = $f[4];
        my $cigar = $f[5];
        my $rnext = $f[6];
        my $pnext = $f[7];
        my $tlen  = $f[8];

        next unless is_paired($flag);
        next if is_unmapped($flag);
        next if is_mate_unmapped($flag);
        next if is_secondary($flag);
        next if is_supplementary($flag);

        if ($exclude_dup) {
            next if is_duplicate($flag);
        }

        next unless defined $mapq && $mapq =~ /^\d+$/;
        next if $mapq < $min_mapq;

        next unless defined $pnext && $pnext =~ /^\d+$/;
        next unless defined $tlen && $tlen =~ /^-?\d+$/;

        my $mate_chr = $rnext eq "=" ? $rname : $rnext;

        # This workflow focuses on intrachromosomal deletion-like events.
        next unless $mate_chr eq $rname;

        my $abs_tlen = abs($tlen);
        next if $abs_tlen < $min_insert_size;

        # Avoid counting both ends of the same read pair.
        next if $tlen < 0;

        my ($left_pos, $right_pos) = sort { $a <=> $b } ($pos, $pnext);

        my $pair_key = join("|",
            $qname,
            $rname,
            $left_pos,
            $right_pos,
            $abs_tlen
        );

        my $pair = {
            sample      => $sample,
            read_name   => $qname,
            chr         => $rname,
            read_pos    => $pos,
            mate_chr    => $mate_chr,
            mate_pos    => $pnext,
            left_pos    => $left_pos,
            right_pos   => $right_pos,
            insert_size => $abs_tlen,
            mapq        => $mapq,
            flag        => $flag,
            cigar       => $cigar,
            target_name => $region->{name},
            pair_key    => $pair_key,
        };

        push @pairs, $pair;

        print $raw_fh join("\t",
            $sample,
            $qname,
            $rname,
            $pos,
            $mate_chr,
            $pnext,
            $abs_tlen,
            $mapq,
            $flag,
            $cigar,
            $region->{name},
            $pair_key,
        ) . "\n";
    }

    close $pipe;

    return @pairs;
}


sub cluster_discordant_pairs {
    my %args = @_;

    my $pairs_ref        = $args{pairs};
    my $cluster_distance = $args{cluster_distance};

    my @pairs = @$pairs_ref;

    my %by_chr;

    foreach my $p (@pairs) {
        push @{ $by_chr{ $p->{chr} } }, $p;
    }

    my @clusters;

    foreach my $chr (sort keys %by_chr) {
        my @sorted = sort {
            $a->{left_pos} <=> $b->{left_pos}
            ||
            $a->{right_pos} <=> $b->{right_pos}
        } @{ $by_chr{$chr} };

        my @current;

        foreach my $p (@sorted) {
            if (!@current) {
                push @current, $p;
                next;
            }

            my $last = $current[-1];

            my $left_close  = abs($p->{left_pos}  - $last->{left_pos})  <= $cluster_distance;
            my $right_close = abs($p->{right_pos} - $last->{right_pos}) <= $cluster_distance;

            if ($left_close && $right_close) {
                push @current, $p;
            }
            else {
                push @clusters, build_cluster(\@current);
                @current = ($p);
            }
        }

        if (@current) {
            push @clusters, build_cluster(\@current);
        }
    }

    return @clusters;
}


sub build_cluster {
    my ($pairs_ref) = @_;

    my @pairs = @$pairs_ref;

    my $first = $pairs[0];

    my @lefts;
    my @rights;
    my @sizes;
    my %target_names;

    foreach my $p (@pairs) {
        push @lefts,  $p->{left_pos};
        push @rights, $p->{right_pos};
        push @sizes,  $p->{insert_size};
        $target_names{ $p->{target_name} } = 1;
    }

    my $left_min  = min(@lefts);
    my $left_max  = max(@lefts);
    my $right_min = min(@rights);
    my $right_max = max(@rights);

    my $candidate_start = $left_max;
    my $candidate_end   = $right_min;

    if ($candidate_end < $candidate_start) {
        $candidate_start = $left_min;
        $candidate_end   = $right_max;
    }

    my $candidate_length = $candidate_end - $candidate_start + 1;
    $candidate_length = 0 if $candidate_length < 0;

    my $target_name = join(",", sort keys %target_names);

    return {
        chr                => $first->{chr},
        start              => $candidate_start,
        end                => $candidate_end,
        candidate_length   => $candidate_length,
        discordant_reads   => scalar(@pairs),
        median_insert_size => median(@sizes),
        min_insert_size    => min(@sizes),
        max_insert_size    => max(@sizes),
        left_min           => $left_min,
        left_max           => $left_max,
        right_min          => $right_min,
        right_max          => $right_max,
        target_name        => $target_name,
        pairs              => \@pairs,
    };
}


sub remove_duplicate_pairs {
    my @pairs = @_;

    my %seen;
    my @unique;

    foreach my $p (@pairs) {
        next if $seen{ $p->{pair_key} }++;
        push @unique, $p;
    }

    return @unique;
}


sub read_target_bed {
    my ($file) = @_;

    my @regions;

    open my $fh, "<", $file
        or die "[ERROR] Cannot open target BED: $file\n";

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        my @f = split /\t/, $line;
        die "[ERROR] Invalid BED line: $line\n" unless @f >= 3;

        my $chr    = $f[0];
        my $start0 = $f[1];
        my $end1   = $f[2];
        my $name   = $f[3] || "$chr:$start0-$end1";

        die "[ERROR] Invalid BED coordinate: $line\n"
            unless $start0 =~ /^\d+$/ && $end1 =~ /^\d+$/ && $end1 > $start0;

        my $start1 = $start0 + 1;

        push @regions, {
            chr    => $chr,
            start  => $start1,
            end    => $end1,
            name   => $name,
            region => "$chr:$start1-$end1",
        };
    }

    close $fh;

    return @regions;
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


sub is_paired {
    my ($flag) = @_;
    return ($flag & 0x1) ? 1 : 0;
}


sub is_unmapped {
    my ($flag) = @_;
    return ($flag & 0x4) ? 1 : 0;
}


sub is_mate_unmapped {
    my ($flag) = @_;
    return ($flag & 0x8) ? 1 : 0;
}


sub is_secondary {
    my ($flag) = @_;
    return ($flag & 0x100) ? 1 : 0;
}


sub is_duplicate {
    my ($flag) = @_;
    return ($flag & 0x400) ? 1 : 0;
}


sub is_supplementary {
    my ($flag) = @_;
    return ($flag & 0x800) ? 1 : 0;
}


sub median {
    my @x = sort { $a <=> $b } @_;

    return 0 unless @x;

    my $n = scalar @x;

    if ($n % 2) {
        return $x[int($n / 2)];
    }
    else {
        return ($x[$n / 2 - 1] + $x[$n / 2]) / 2;
    }
}


sub min {
    my @x = @_;
    return 0 unless @x;

    my $min = $x[0];

    foreach my $v (@x) {
        $min = $v if $v < $min;
    }

    return $min;
}


sub max {
    my @x = @_;
    return 0 unless @x;

    my $max = $x[0];

    foreach my $v (@x) {
        $max = $v if $v > $max;
    }

    return $max;
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

  perl bin/run_discordant_reads.pl \\
    --config conf/hcm_exondel.conf \\
    --bam sample.sorted.bam \\
    --sample SAMPLE001 \\
    --out results/SAMPLE001/02.discordant_reads/SAMPLE001.discordant_reads.tsv

Required command-line arguments:

  --config       Config file
  --bam          Coordinate-sorted BAM file
  --sample       Sample ID
  --out          Output discordant read summary file

All analysis parameters are read from config file:

  SAMTOOLS
  TARGET_REGION_BED
  SCAN_WHOLE_BAM
  MIN_MAPQ
  MIN_DISCORDANT_INSERT_SIZE
  MIN_DISCORDANT_READS
  DISCORDANT_CLUSTER_DISTANCE
  MIN_DISCORDANT_CANDIDATE_SIZE
  MAX_DISCORDANT_CANDIDATE_SIZE
  EXCLUDE_DUPLICATES

Default behavior:

  SCAN_WHOLE_BAM=0

  The script scans only candidate gene regions defined by TARGET_REGION_BED.

Output:

  *.discordant_reads.tsv
  *.discordant_reads.raw_discordant_pairs.tsv
  *.discordant_reads.supporting_pairs.tsv

USAGE
}

