#!/usr/bin/env perl
# -*- coding: utf-8 -*-

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path);

# ============================================================
# run_split_reads.pl
#
# Split-read / soft-clipped read screening from coordinate-sorted
# WGS BAM.
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
#   MIN_SOFTCLIP_LENGTH
#   MIN_SPLIT_READS
#   SPLIT_CLUSTER_DISTANCE
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

my $scan_whole_bam = normalize_bool(
    get_conf_value(\%CONF, "SCAN_WHOLE_BAM", 0)
);

my $target_bed = get_conf_value(\%CONF, "TARGET_REGION_BED", "");

if (!$scan_whole_bam) {
    die "[ERROR] TARGET_REGION_BED is required in config when SCAN_WHOLE_BAM=0\n"
        unless $target_bed;

    $target_bed = abs_path($target_bed);
    die "[ERROR] TARGET_REGION_BED file not found: $target_bed\n"
        unless -s $target_bed;
}

my $min_mapq = get_conf_value(\%CONF, "MIN_MAPQ", 20);

my $min_softclip_length =
    get_conf_value(\%CONF, "MIN_SOFTCLIP_LENGTH", 10);

my $min_split_reads =
    get_conf_value(\%CONF, "MIN_SPLIT_READS", 2);

my $split_cluster_distance =
    get_conf_value(\%CONF, "SPLIT_CLUSTER_DISTANCE", 200);

my $exclude_duplicates =
    normalize_bool(get_conf_value(\%CONF, "EXCLUDE_DUPLICATES", 1));

# ============================================================
# Prepare output
# ============================================================

my $outdir = dirname($out);
make_path($outdir) unless -d $outdir;

my $prefix = $out;
$prefix =~ s/\.tsv$//;

my $raw_split_reads_out = $prefix . ".raw_split_reads.tsv";
my $supporting_reads_out = $prefix . ".supporting_reads.tsv";

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

print "[INFO] Split-read screening started\n";
print "[INFO] Sample                  : $sample\n";
print "[INFO] BAM                     : $bam\n";
print "[INFO] Config                  : $config\n";
print "[INFO] Scan whole BAM           : $scan_whole_bam\n";
print "[INFO] Target BED               : " . ($target_bed || "NA") . "\n";
print "[INFO] Scan region number       : " . scalar(@scan_regions) . "\n";
print "[INFO] Min MAPQ                 : $min_mapq\n";
print "[INFO] Min soft-clip length     : $min_softclip_length\n";
print "[INFO] Min split reads          : $min_split_reads\n";
print "[INFO] Split cluster distance   : $split_cluster_distance\n";

# ============================================================
# Scan BAM for split reads
# ============================================================

my @split_reads;

open my $raw_fh, ">", $raw_split_reads_out
    or die "[ERROR] Cannot write raw split reads output: $raw_split_reads_out\n";

print $raw_fh join("\t",
    qw(
        SampleID
        Read_Name
        Chrom
        Pos
        Breakpoint_Pos
        MAPQ
        FLAG
        CIGAR
        Softclip_Side
        Softclip_Length
        Has_SA_Tag
        SA_Tag
        Target_Name
        Read_Key
    )
) . "\n";

foreach my $region (@scan_regions) {
    my @reads = scan_region_for_split_reads(
        bam                 => $bam,
        samtools            => $samtools,
        sample              => $sample,
        region              => $region,
        min_mapq            => $min_mapq,
        min_softclip_length => $min_softclip_length,
        exclude_dup         => $exclude_duplicates,
        raw_fh              => $raw_fh,
    );

    push @split_reads, @reads;
}

close $raw_fh;

my @unique_split_reads = remove_duplicate_reads(@split_reads);

# ============================================================
# Cluster split reads
# ============================================================

my @clusters = cluster_split_reads(
    reads            => \@unique_split_reads,
    cluster_distance => $split_cluster_distance,
);

# ============================================================
# Output summary and supporting reads
# ============================================================

open my $out_fh, ">", $out
    or die "[ERROR] Cannot write output file: $out\n";

open my $support_fh, ">", $supporting_reads_out
    or die "[ERROR] Cannot write supporting reads output: $supporting_reads_out\n";

print $out_fh join("\t",
    qw(
        SampleID
        Chrom
        Start
        End
        Candidate_Length
        Split_Reads
        Left_Softclip_Reads
        Right_Softclip_Reads
        SA_Tag_Reads
        Breakpoint_Min
        Breakpoint_Max
        Target_Name
        Split_Status
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
        Pos
        Breakpoint_Pos
        MAPQ
        FLAG
        CIGAR
        Softclip_Side
        Softclip_Length
        Has_SA_Tag
        SA_Tag
        Target_Name
    )
) . "\n";

my $cluster_id = 0;
my $pass_cluster_count = 0;

foreach my $cluster (@clusters) {
    next if $cluster->{split_reads} < $min_split_reads;

    $cluster_id++;
    $pass_cluster_count++;

    my $status = "Split_supported";
    my $candidate_status = "Split_candidate";
    my $comment = "Candidate breakpoint region supported by split or soft-clipped reads";

    print $out_fh join("\t",
        $sample,
        $cluster->{chr},
        $cluster->{start},
        $cluster->{end},
        $cluster->{candidate_length},
        $cluster->{split_reads},
        $cluster->{left_softclip_reads},
        $cluster->{right_softclip_reads},
        $cluster->{sa_tag_reads},
        $cluster->{breakpoint_min},
        $cluster->{breakpoint_max},
        $cluster->{target_name},
        $status,
        $candidate_status,
        $comment,
    ) . "\n";

    foreach my $r (@{ $cluster->{reads} }) {
        print $support_fh join("\t",
            $sample,
            "SR_CLUSTER_" . $cluster_id,
            $r->{read_name},
            $r->{chr},
            $r->{pos},
            $r->{breakpoint_pos},
            $r->{mapq},
            $r->{flag},
            $r->{cigar},
            $r->{softclip_side},
            $r->{softclip_length},
            $r->{has_sa_tag},
            $r->{sa_tag},
            $r->{target_name},
        ) . "\n";
    }
}

close $out_fh;
close $support_fh;

print "[INFO] Split-read screening finished\n";
print "[INFO] Raw split reads output    : $raw_split_reads_out\n";
print "[INFO] Summary output            : $out\n";
print "[INFO] Supporting reads output   : $supporting_reads_out\n";
print "[INFO] Raw split reads           : " . scalar(@split_reads) . "\n";
print "[INFO] Unique split reads        : " . scalar(@unique_split_reads) . "\n";
print "[INFO] Total clusters            : " . scalar(@clusters) . "\n";
print "[INFO] Passed clusters           : $pass_cluster_count\n";

exit 0;


# ============================================================
# Subroutines
# ============================================================

sub scan_region_for_split_reads {
    my %args = @_;

    my $bam                 = $args{bam};
    my $samtools            = $args{samtools};
    my $sample              = $args{sample};
    my $region              = $args{region};
    my $min_mapq            = $args{min_mapq};
    my $min_softclip_length = $args{min_softclip_length};
    my $exclude_dup         = $args{exclude_dup};
    my $raw_fh              = $args{raw_fh};

    my @reads;

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

        next if is_unmapped($flag);
        next if is_secondary($flag);

        if ($exclude_dup) {
            next if is_duplicate($flag);
        }

        next unless defined $mapq && $mapq =~ /^\d+$/;
        next if $mapq < $min_mapq;

        next if !$cigar || $cigar eq "*";

        my ($softclip_side, $softclip_len, $breakpoint_pos) =
            parse_softclip_breakpoint($pos, $cigar);

        my ($has_sa_tag, $sa_tag) = parse_sa_tag(@f);

        my $is_supplementary = is_supplementary($flag);

        # Keep reads with sufficient soft clipping, SA tag, or supplementary flag.
        my $keep = 0;

        if ($softclip_len >= $min_softclip_length) {
            $keep = 1;
        }

        if ($has_sa_tag) {
            $keep = 1;
        }

        if ($is_supplementary) {
            $keep = 1;
        }

        next unless $keep;

        if ($softclip_len < $min_softclip_length && !$has_sa_tag && !$is_supplementary) {
            next;
        }

        if (!$breakpoint_pos || $breakpoint_pos < 1) {
            $breakpoint_pos = $pos;
        }

        my $read_key = join("|",
            $qname,
            $rname,
            $pos,
            $cigar,
            $breakpoint_pos
        );

        my $read = {
            sample          => $sample,
            read_name       => $qname,
            chr             => $rname,
            pos             => $pos,
            breakpoint_pos  => $breakpoint_pos,
            mapq            => $mapq,
            flag            => $flag,
            cigar           => $cigar,
            softclip_side   => $softclip_side,
            softclip_length => $softclip_len,
            has_sa_tag      => $has_sa_tag,
            sa_tag          => $sa_tag,
            target_name     => $region->{name},
            read_key        => $read_key,
        };

        push @reads, $read;

        print $raw_fh join("\t",
            $sample,
            $qname,
            $rname,
            $pos,
            $breakpoint_pos,
            $mapq,
            $flag,
            $cigar,
            $softclip_side,
            $softclip_len,
            $has_sa_tag,
            $sa_tag,
            $region->{name},
            $read_key,
        ) . "\n";
    }

    close $pipe;

    return @reads;
}


sub parse_softclip_breakpoint {
    my ($pos, $cigar) = @_;

    my $left_softclip  = 0;
    my $right_softclip = 0;

    if ($cigar =~ /^(\d+)S/) {
        $left_softclip = $1;
    }

    if ($cigar =~ /(\d+)S$/) {
        $right_softclip = $1;
    }

    my $ref_len = cigar_ref_length($cigar);
    my $align_end = $pos + $ref_len - 1;

    if ($left_softclip >= $right_softclip && $left_softclip > 0) {
        return ("Left", $left_softclip, $pos);
    }
    elsif ($right_softclip > 0) {
        return ("Right", $right_softclip, $align_end);
    }
    else {
        return ("None", 0, $pos);
    }
}


sub parse_sa_tag {
    my (@fields) = @_;

    my $sa_tag = "NA";

    for (my $i = 11; $i < @fields; $i++) {
        if ($fields[$i] =~ /^SA:Z:(.+)$/) {
            $sa_tag = $1;
            return (1, $sa_tag);
        }
    }

    return (0, $sa_tag);
}


sub cigar_ref_length {
    my ($cigar) = @_;

    my $ref_len = 0;

    while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
        my $len = $1;
        my $op  = $2;

        if ($op =~ /[MDN=X]/) {
            $ref_len += $len;
        }
    }

    return $ref_len;
}


sub cluster_split_reads {
    my %args = @_;

    my $reads_ref        = $args{reads};
    my $cluster_distance = $args{cluster_distance};

    my @reads = @$reads_ref;

    my %by_chr;

    foreach my $r (@reads) {
        push @{ $by_chr{ $r->{chr} } }, $r;
    }

    my @clusters;

    foreach my $chr (sort keys %by_chr) {
        my @sorted = sort {
            $a->{breakpoint_pos} <=> $b->{breakpoint_pos}
        } @{ $by_chr{$chr} };

        my @current;

        foreach my $r (@sorted) {
            if (!@current) {
                push @current, $r;
                next;
            }

            my $last = $current[-1];

            if (abs($r->{breakpoint_pos} - $last->{breakpoint_pos}) <= $cluster_distance) {
                push @current, $r;
            }
            else {
                push @clusters, build_split_cluster(\@current);
                @current = ($r);
            }
        }

        if (@current) {
            push @clusters, build_split_cluster(\@current);
        }
    }

    return @clusters;
}


sub build_split_cluster {
    my ($reads_ref) = @_;

    my @reads = @$reads_ref;

    my $first = $reads[0];

    my @breakpoints;
    my %target_names;

    my $left_softclip_reads  = 0;
    my $right_softclip_reads = 0;
    my $sa_tag_reads         = 0;

    foreach my $r (@reads) {
        push @breakpoints, $r->{breakpoint_pos};
        $target_names{ $r->{target_name} } = 1;

        if ($r->{softclip_side} eq "Left") {
            $left_softclip_reads++;
        }
        elsif ($r->{softclip_side} eq "Right") {
            $right_softclip_reads++;
        }

        if ($r->{has_sa_tag}) {
            $sa_tag_reads++;
        }
    }

    my $bp_min = min(@breakpoints);
    my $bp_max = max(@breakpoints);

    my $candidate_length = $bp_max - $bp_min + 1;
    $candidate_length = 0 if $candidate_length < 0;

    my $target_name = join(",", sort keys %target_names);

    return {
        chr                  => $first->{chr},
        start                => $bp_min,
        end                  => $bp_max,
        candidate_length     => $candidate_length,
        split_reads          => scalar(@reads),
        left_softclip_reads  => $left_softclip_reads,
        right_softclip_reads => $right_softclip_reads,
        sa_tag_reads         => $sa_tag_reads,
        breakpoint_min       => $bp_min,
        breakpoint_max       => $bp_max,
        target_name          => $target_name,
        reads                => \@reads,
    };
}


sub remove_duplicate_reads {
    my @reads = @_;

    my %seen;
    my @unique;

    foreach my $r (@reads) {
        next if $seen{ $r->{read_key} }++;
        push @unique, $r;
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


sub is_unmapped {
    my ($flag) = @_;
    return ($flag & 0x4) ? 1 : 0;
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

  perl bin/run_split_reads.pl \\
    --config conf/hcm_exondel.conf \\
    --bam sample.sorted.bam \\
    --sample SAMPLE001 \\
    --out results/SAMPLE001/03.split_reads/SAMPLE001.split_reads.tsv

Required command-line arguments:

  --config       Config file
  --bam          Coordinate-sorted BAM file
  --sample       Sample ID
  --out          Output split-read summary file

All analysis parameters are read from config file:

  SAMTOOLS
  TARGET_REGION_BED
  SCAN_WHOLE_BAM
  MIN_MAPQ
  MIN_SOFTCLIP_LENGTH
  MIN_SPLIT_READS
  SPLIT_CLUSTER_DISTANCE
  EXCLUDE_DUPLICATES

Default behavior:

  SCAN_WHOLE_BAM=0

  The script scans only candidate gene regions defined by TARGET_REGION_BED.

Output:

  *.split_reads.tsv
  *.split_reads.raw_split_reads.tsv
  *.split_reads.supporting_reads.tsv

USAGE
}

