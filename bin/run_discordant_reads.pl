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
# Optimized according to conf/hcm_exondel.example.conf:
#   1. Target regions are read from REFSEQ_MANE_SELECT_GENE_TXT.
#   2. All TXT interval coordinates are 1-based closed.
#   3. No BED file is required.
#   4. TARGET_REGION_FLANK is applied internally when SCAN_WHOLE_BAM=0.
#   5. MIN_DISCORDANT_MAPQ is used for discordant read-pair filtering,
#      with fallback to MIN_MAPQ.
#   6. Cluster filtering only uses MIN_DISCORDANT_READS.
#      MIN_DISCORDANT_CANDIDATE_SIZE and MAX_DISCORDANT_CANDIDATE_SIZE
#      are intentionally not used.
#   7. Relative paths in config such as db/... are resolved against
#      the project root inferred from the config path.
#
# Required command-line arguments:
#   --config
#   --bam
#   --sample
#   --out
#
# Main config-controlled parameters:
#   SAMTOOLS
#   REFSEQ_MANE_SELECT_GENE_TXT
#   REF_FASTA_INDEX
#   SCAN_WHOLE_BAM
#   TARGET_REGION_FLANK
#   MIN_MAPQ
#   MIN_DISCORDANT_MAPQ
#   MIN_DISCORDANT_INSERT_SIZE
#   MIN_DISCORDANT_READS
#   DISCORDANT_CLUSTER_DISTANCE
#   EXCLUDE_DUPLICATES
#   VERBOSE
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

die "[ERROR] Config file not found: $config\n" unless defined $config && -s $config;
die "[ERROR] BAM file not found: $bam\n" unless defined $bam && -s $bam;

my %CONF = read_config($config);

# Project root inference:
# If config is /path/to/HCM_WGS_EXON_DEL/conf/hcm_exondel.conf,
# project_root is /path/to/HCM_WGS_EXON_DEL.
# Therefore config values like db/xxx are resolved as project_root/db/xxx.
my $config_dir   = dirname($config);
my $project_root = dirname($config_dir);

# ============================================================
# Read parameters from config
# ============================================================
my $samtools = get_conf_required(\%CONF, "SAMTOOLS");

# Keep command names such as "samtools" unchanged so that PATH can resolve them.
# Only absolute SAMTOOLS paths are checked here.
if ($samtools =~ m{^/}) {
    die "[ERROR] SAMTOOLS not executable: $samtools\n" unless -x $samtools;
}

my $scan_whole_bam = normalize_bool(get_conf_value(\%CONF, "SCAN_WHOLE_BAM", 0));
my $target_flank   = get_conf_value(\%CONF, "TARGET_REGION_FLANK", 0);
validate_nonnegative_integer("TARGET_REGION_FLANK", $target_flank);

my $min_mapq = get_conf_value(
    \%CONF,
    "MIN_DISCORDANT_MAPQ",
    get_conf_value(\%CONF, "MIN_MAPQ", 20)
);

my $min_discordant_insert_size = get_conf_value(\%CONF, "MIN_DISCORDANT_INSERT_SIZE", 1000);
my $min_discordant_reads       = get_conf_value(\%CONF, "MIN_DISCORDANT_READS", 3);
my $cluster_distance           = get_conf_value(\%CONF, "DISCORDANT_CLUSTER_DISTANCE", 1000);
my $exclude_duplicates         = normalize_bool(get_conf_value(\%CONF, "EXCLUDE_DUPLICATES", 1));

# For deletion-like discordant read pairs, only keep inward-facing FR orientation:
#   left read  = Forward
#   right read = Reverse
# This keeps both F1R2 and F2R1, and filters R1F2/R2F1.
# Default: 1
my $filter_deletion_orientation = normalize_bool(get_conf_value(\%CONF, "FILTER_DELETION_ORIENTATION", 1));

my $verbose = normalize_bool(get_conf_value(\%CONF, "VERBOSE", 1));

validate_nonnegative_integer("MIN_DISCORDANT_MAPQ/MIN_MAPQ", $min_mapq);
validate_positive_integer("MIN_DISCORDANT_INSERT_SIZE", $min_discordant_insert_size);
validate_positive_integer("MIN_DISCORDANT_READS", $min_discordant_reads);
validate_positive_integer("DISCORDANT_CLUSTER_DISTANCE", $cluster_distance);

my $gene_txt = "";
my $hcm_gene_list = "";
my %target_gene_set;
my %chr_length;

if (!$scan_whole_bam) {
    $gene_txt = get_conf_required(\%CONF, "REFSEQ_MANE_SELECT_GENE_TXT");
    $gene_txt = resolve_config_path($gene_txt, $project_root);

    die "[ERROR] REFSEQ_MANE_SELECT_GENE_TXT file not found: $gene_txt
" unless -s $gene_txt;

    # When SCAN_WHOLE_BAM=0, restrict target scanning to genes in HCM_CORE_GENE_LIST if configured.
    # REFSEQ_MANE_SELECT_GENE_TXT remains the coordinate source, while HCM_CORE_GENE_LIST defines
    # which genes should be scanned.
    $hcm_gene_list = get_conf_value(\%CONF, "HCM_CORE_GENE_LIST", "");
    if ($hcm_gene_list) {
        $hcm_gene_list = resolve_config_path($hcm_gene_list, $project_root);
        die "[ERROR] HCM_CORE_GENE_LIST file not found: $hcm_gene_list
" unless -s $hcm_gene_list;
        %target_gene_set = read_gene_list($hcm_gene_list);
        die "[ERROR] No valid genes found in HCM_CORE_GENE_LIST: $hcm_gene_list
" unless scalar(keys %target_gene_set) > 0;
    }

    my $fai = get_conf_value(\%CONF, "REF_FASTA_INDEX", "");
    if ($fai) {
        $fai = resolve_config_path($fai, $project_root);
        if (-s $fai) {
            %chr_length = read_fai($fai);
        } else {
            warn "[WARN] REF_FASTA_INDEX is configured but not found: $fai. Region ends will not be clipped by chromosome length.\n";
        }
    }
}

# ============================================================
# Prepare output
# ============================================================
my $outdir = dirname($out);
make_path($outdir) unless -d $outdir;

my $prefix = $out;
$prefix =~ s/\.tsv$//;

my $supporting_pairs_out = $prefix . ".supporting_pairs.tsv";
my $raw_pairs_out        = $prefix . ".raw_discordant_pairs.tsv";
my $discard_log_out      = $prefix . ".discarded_clusters.tsv";

# ============================================================
# Prepare scan regions
# ============================================================
my @scan_regions;

if ($scan_whole_bam) {
    push @scan_regions, {
        gene       => "ALL",
        transcript => "NA",
        chr        => "ALL",
        start      => 0,
        end        => 0,
        strand     => "NA",
        exon_count => "NA",
        scan_start => 0,
        scan_end   => 0,
        name       => "whole_bam",
        region     => "",
    };
} else {
    @scan_regions = read_gene_txt(
        file            => $gene_txt,
        flank           => $target_flank,
        chr_length      => \%chr_length,
        target_gene_set => \%target_gene_set,
    );
    die "[ERROR] No valid target regions found in REFSEQ_MANE_SELECT_GENE_TXT: $gene_txt\n" unless @scan_regions;
}

log_msg($verbose, "[INFO] Discordant read screening started");
log_msg($verbose, "[INFO] Sample: $sample");
log_msg($verbose, "[INFO] BAM: $bam");
log_msg($verbose, "[INFO] Config: $config");
log_msg($verbose, "[INFO] Project root: $project_root");
log_msg($verbose, "[INFO] Scan whole BAM: $scan_whole_bam");
log_msg($verbose, "[INFO] Gene TXT: " . ($gene_txt || "NA"));
log_msg($verbose, "[INFO] HCM core gene list: " . ($hcm_gene_list || "NA"));
log_msg($verbose, "[INFO] HCM core gene number: " . scalar(keys %target_gene_set));
log_msg($verbose, "[INFO] Target flank: $target_flank");
log_msg($verbose, "[INFO] Scan region number: " . scalar(@scan_regions));
log_msg($verbose, "[INFO] Min discordant MAPQ: $min_mapq");
log_msg($verbose, "[INFO] Min discordant insert size: $min_discordant_insert_size");
log_msg($verbose, "[INFO] Min discordant reads: $min_discordant_reads");
log_msg($verbose, "[INFO] Cluster distance: $cluster_distance");
log_msg($verbose, "[INFO] Exclude duplicates: $exclude_duplicates");
log_msg($verbose, "[INFO] Filter deletion-like FR orientation: $filter_deletion_orientation");

# ============================================================
# Scan BAM for discordant pairs
# ============================================================
my @discordant_pairs;

open my $raw_fh, ">", $raw_pairs_out
    or die "[ERROR] Cannot write raw pairs output: $raw_pairs_out\n";

print $raw_fh join("    ", qw(
    SampleID Read_Name Chrom Read_Pos Mate_Chrom Mate_Pos Insert_Size
    MAPQ FLAG CIGAR Pair_Orientation Gene Transcript Target_Name Scan_Region Pair_Key
)) . "
";

foreach my $region (@scan_regions) {
    my @pairs = scan_region_for_discordant_pairs(
        bam             => $bam,
        samtools        => $samtools,
        sample          => $sample,
        region          => $region,
        min_mapq        => $min_mapq,
        min_insert_size => $min_discordant_insert_size,
        exclude_dup     => $exclude_duplicates,
        filter_del_ori => $filter_deletion_orientation,
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
# Output summary, supporting reads, and discarded clusters
# ============================================================
open my $out_fh, ">", $out
    or die "[ERROR] Cannot write output file: $out\n";

open my $support_fh, ">", $supporting_pairs_out
    or die "[ERROR] Cannot write supporting pairs output: $supporting_pairs_out\n";

open my $discard_fh, ">", $discard_log_out
    or die "[ERROR] Cannot write discarded cluster log: $discard_log_out\n";

print $out_fh join("\t", qw(
    SampleID Cluster_ID Chrom Start End Candidate_Length Discordant_Reads
    Median_Insert_Size Min_Insert_Size Max_Insert_Size
    Left_Boundary_Min Left_Boundary_Max Right_Boundary_Min Right_Boundary_Max
    Gene Transcript Target_Name Discordant_Status Candidate_Status Comment
)) . "\n";

print $support_fh join("    ", qw(
    SampleID Cluster_ID Read_Name Chrom Read_Pos Mate_Chrom Mate_Pos Insert_Size
    MAPQ FLAG CIGAR Pair_Orientation Gene Transcript Target_Name Scan_Region
)) . "
";

print $discard_fh join("\t", qw(
    SampleID Chrom Start End Candidate_Length Discordant_Reads
    Median_Insert_Size Gene Transcript Target_Name Discard_Reason
)) . "\n";

my $cluster_id = 0;
my $pass_cluster_count = 0;
my $discard_cluster_count = 0;

foreach my $cluster (@clusters) {
    my @discard_reasons;
    push @discard_reasons, "low_discordant_reads" if $cluster->{discordant_reads} < $min_discordant_reads;

    if (@discard_reasons) {
        $discard_cluster_count++;
        print $discard_fh join("\t",
            $sample,
            $cluster->{chr},
            $cluster->{start},
            $cluster->{end},
            $cluster->{candidate_length},
            $cluster->{discordant_reads},
            sprintf("%.2f", $cluster->{median_insert_size}),
            $cluster->{gene},
            $cluster->{transcript},
            $cluster->{target_name},
            join(",", @discard_reasons),
        ) . "\n";
        next;
    }

    $cluster_id++;
    $pass_cluster_count++;

    my $cluster_name      = "DR_CLUSTER_" . $cluster_id;
    my $status            = "Discordant_supported";
    my $candidate_status  = "Discordant_candidate";
    my $comment           = "Deletion-like candidate supported by intrachromosomal discordant read pairs";

    print $out_fh join("\t",
        $sample,
        $cluster_name,
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
        $cluster->{gene},
        $cluster->{transcript},
        $cluster->{target_name},
        $status,
        $candidate_status,
        $comment,
    ) . "\n";

    foreach my $p (@{ $cluster->{pairs} }) {
        print $support_fh join("\t",
            $sample,
            $cluster_name,
            $p->{read_name},
            $p->{chr},
            $p->{read_pos},
            $p->{mate_chr},
            $p->{mate_pos},
            $p->{insert_size},
            $p->{mapq},
            $p->{flag},
            $p->{cigar},
            $p->{pair_orientation},
            $p->{gene},
            $p->{transcript},
            $p->{target_name},
            $p->{scan_region},
        ) . "\n";
    }
}

close $out_fh;
close $support_fh;
close $discard_fh;

log_msg($verbose, "[INFO] Discordant reads screening finished");
log_msg($verbose, "[INFO] Raw discordant pairs output: $raw_pairs_out");
log_msg($verbose, "[INFO] Summary output: $out");
log_msg($verbose, "[INFO] Supporting pairs output: $supporting_pairs_out");
log_msg($verbose, "[INFO] Discarded cluster log: $discard_log_out");
log_msg($verbose, "[INFO] Raw discordant pairs: " . scalar(@discordant_pairs));
log_msg($verbose, "[INFO] Unique discordant pairs: " . scalar(@unique_pairs));
log_msg($verbose, "[INFO] Total clusters: " . scalar(@clusters));
log_msg($verbose, "[INFO] Passed clusters: $pass_cluster_count");
log_msg($verbose, "[INFO] Discarded clusters: $discard_cluster_count");

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
    my $filter_del_ori = $args{filter_del_ori};
    my $raw_fh          = $args{raw_fh};

    my @pairs;
    my $cmd;

    if ($region->{region}) {
        $cmd = join(" ", shell_quote($samtools), "view", shell_quote($bam), shell_quote($region->{region}));
    } else {
        $cmd = join(" ", shell_quote($samtools), "view", shell_quote($bam));
    }

    open my $pipe, "$cmd |" or die "[ERROR] Failed to run command: $cmd\n";

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

        next unless defined $flag && $flag =~ /^\d+$/;
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

        next unless defined $pos   && $pos   =~ /^\d+$/;
        next unless defined $pnext && $pnext =~ /^\d+$/;
        next unless defined $tlen  && $tlen  =~ /^-?\d+$/;

        my $mate_chr = $rnext eq "=" ? $rname : $rnext;

        # This module focuses on intrachromosomal deletion-like events.
        next unless $mate_chr eq $rname;

        my $pair_orientation = get_pair_orientation($flag, $pos, $pnext);

        # Deletion-like discordant read pairs should have FR orientation when ordered by genomic position:
        #   left read  = Forward
        #   right read = Reverse
        # This keeps both F1R2 and F2R1.
        if ($filter_del_ori) {
            next unless is_deletion_like_orientation($flag, $pos, $pnext);
        }

        my $abs_tlen = abs($tlen);
        next if $abs_tlen < $min_insert_size;

        # Keep only one side of the read pair to avoid double-counting.
        next if $tlen < 0;

        my ($left_pos, $right_pos) = sort { $a <=> $b } ($pos, $pnext);

        my $pair_key = join("|",
            $qname,
            $rname,
            $left_pos,
            $right_pos,
            $abs_tlen,
        );

        my $scan_region = $region->{region} || "whole_bam";

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
            pair_orientation => $pair_orientation,
            gene        => $region->{gene},
            transcript  => $region->{transcript},
            target_name => $region->{name},
            scan_region => $scan_region,
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
            $pair_orientation,
            $region->{gene},
            $region->{transcript},
            $region->{name},
            $scan_region,
            $pair_key,
        ) . "\n";
    }

    my $ok = close $pipe;
    die "[ERROR] samtools view command failed: $cmd\n" unless $ok;

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
            $a->{left_pos}  <=> $b->{left_pos} ||
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
            } else {
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
    my %genes;
    my %transcripts;
    my %target_names;

    foreach my $p (@pairs) {
        push @lefts,  $p->{left_pos};
        push @rights, $p->{right_pos};
        push @sizes,  $p->{insert_size};
        $genes{ $p->{gene} }               = 1 if defined $p->{gene};
        $transcripts{ $p->{transcript} }   = 1 if defined $p->{transcript};
        $target_names{ $p->{target_name} } = 1 if defined $p->{target_name};
    }

    my $left_min  = min(@lefts);
    my $left_max  = max(@lefts);
    my $right_min = min(@rights);
    my $right_max = max(@rights);

    # Conservative interval between inner boundaries.
    my $candidate_start = $left_max;
    my $candidate_end   = $right_min;

    # Fallback for abnormal or overlapping signals.
    if ($candidate_end < $candidate_start) {
        $candidate_start = $left_min;
        $candidate_end   = $right_max;
    }

    my $candidate_length = $candidate_end - $candidate_start + 1;
    $candidate_length = 0 if $candidate_length < 0;

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
        gene               => join(",", sort keys %genes),
        transcript         => join(",", sort keys %transcripts),
        target_name        => join(",", sort keys %target_names),
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

sub read_gene_txt {
    my %args = @_;

    my $file           = $args{file};
    my $flank          = $args{flank};
    my $chr_length_ref = $args{chr_length} || {};
    my $target_gene_set = $args{target_gene_set} || {};

    my @regions;

    open my $fh, "<", $file or die "[ERROR] Cannot open gene TXT: $file\n";

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        my @f = split /\t/, $line;

        # Expected columns:
        # Gene Transcript Chrom Start End Strand ExonCount
        die "[ERROR] Invalid gene TXT line. Expected >=7 columns: $line\n" unless @f >= 7;

        my ($gene, $transcript, $chr, $start, $end, $strand, $exon_count) = @f[0..6];

        # Skip a header line if present.
        next if $gene =~ /^Gene$/i && $transcript =~ /^Transcript$/i;

        # If HCM_CORE_GENE_LIST is configured, only keep genes in that list.
        if (scalar(keys %$target_gene_set) > 0) {
            next unless exists $target_gene_set->{$gene};
        }

        die "[ERROR] Invalid gene TXT coordinate: $line\n"
            unless defined $start && defined $end && $start =~ /^\d+$/ && $end =~ /^\d+$/ && $end >= $start;

        my $scan_start = $start - $flank;
        $scan_start = 1 if $scan_start < 1;

        my $scan_end = $end + $flank;
        if (exists $chr_length_ref->{$chr} && $scan_end > $chr_length_ref->{$chr}) {
            $scan_end = $chr_length_ref->{$chr};
        }

        my $name = join("|", $gene, $transcript, $chr . ":" . $start . "-" . $end);
        my $region = $chr . ":" . $scan_start . "-" . $scan_end;

        push @regions, {
            gene       => $gene,
            transcript => $transcript,
            chr        => $chr,
            start      => $start,
            end        => $end,
            strand     => $strand,
            exon_count => $exon_count,
            scan_start => $scan_start,
            scan_end   => $scan_end,
            name       => $name,
            region     => $region,
        };
    }

    close $fh;
    return @regions;
}

sub read_gene_list {
    my ($file) = @_;
    my %genes;

    open my $fh, "<", $file or die "[ERROR] Cannot open gene list: $file
";

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        $line =~ s/^\s+//;
        $line =~ s/\s+$//;
        next if $line eq "";
        next if $line =~ /^#/;

        my @f = split /[\t, ]+/, $line;
        my $gene = $f[0];
        next unless defined $gene && $gene ne "";
        next if $gene =~ /^Gene$/i;

        $genes{$gene} = 1;
    }

    close $fh;
    return %genes;
}

sub read_fai {
    my ($file) = @_;
    my %len;

    open my $fh, "<", $file or die "[ERROR] Cannot open FASTA index: $file\n";

    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^\s*$/;
        my @f = split /\t/, $line;
        next unless @f >= 2;
        my ($chr, $length) = @f[0,1];
        next unless defined $length && $length =~ /^\d+$/;
        $len{$chr} = $length;
    }

    close $fh;
    return %len;
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
            $val =~ s/\s+#.*$//;
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
        unless exists $conf_ref->{$key} && defined $conf_ref->{$key} && $conf_ref->{$key} ne "";
    return $conf_ref->{$key};
}

sub get_conf_value {
    my ($conf_ref, $key, $default) = @_;
    if (exists $conf_ref->{$key} && defined $conf_ref->{$key} && $conf_ref->{$key} ne "") {
        return $conf_ref->{$key};
    }
    return $default;
}

sub resolve_config_path {
    my ($path, $project_root) = @_;

    return "" unless defined $path && $path ne "";

    # Absolute path.
    return $path if $path =~ m{^/};

    # Keep simple command names unchanged. This is mostly for executable fields,
    # but it is safe because file fields are checked by -s later.
    # For paths containing /, resolve them against project root.
    if ($path !~ m{/}) {
        return $path;
    }

    my $resolved = $project_root . "/" . $path;
    $resolved =~ s{//+}{/}g;
    return $resolved;
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

sub validate_positive_integer {
    my ($name, $value) = @_;
    die "[ERROR] $name must be a positive integer. Current value: $value\n"
        unless defined $value && $value =~ /^[1-9]\d*$/;
}

sub validate_nonnegative_integer {
    my ($name, $value) = @_;
    die "[ERROR] $name must be a non-negative integer. Current value: $value\n"
        unless defined $value && $value =~ /^\d+$/;
}

sub is_paired        { my ($flag) = @_; return ($flag & 0x1)   ? 1 : 0; }
sub is_unmapped      { my ($flag) = @_; return ($flag & 0x4)   ? 1 : 0; }
sub is_mate_unmapped { my ($flag) = @_; return ($flag & 0x8)   ? 1 : 0; }
sub is_reverse       { my ($flag) = @_; return ($flag & 0x10)  ? 1 : 0; }
sub is_mate_reverse  { my ($flag) = @_; return ($flag & 0x20)  ? 1 : 0; }
sub is_read1         { my ($flag) = @_; return ($flag & 0x40)  ? 1 : 0; }
sub is_read2         { my ($flag) = @_; return ($flag & 0x80)  ? 1 : 0; }
sub is_secondary     { my ($flag) = @_; return ($flag & 0x100) ? 1 : 0; }
sub is_duplicate     { my ($flag) = @_; return ($flag & 0x400) ? 1 : 0; }
sub is_supplementary { my ($flag) = @_; return ($flag & 0x800) ? 1 : 0; }

sub is_deletion_like_orientation {
    my ($flag, $pos, $pnext) = @_;

    my $read_reverse = is_reverse($flag);
    my $mate_reverse = is_mate_reverse($flag);

    my $left_is_read = ($pos <= $pnext) ? 1 : 0;

    my ($left_reverse, $right_reverse);
    if ($left_is_read) {
        $left_reverse  = $read_reverse;
        $right_reverse = $mate_reverse;
    } else {
        $left_reverse  = $mate_reverse;
        $right_reverse = $read_reverse;
    }

    # FR orientation after sorting by genomic coordinate:
    #   left read  is forward  => left_reverse  == 0
    #   right read is reverse  => right_reverse == 1
    return (!$left_reverse && $right_reverse) ? 1 : 0;
}

sub get_pair_orientation {
    my ($flag, $pos, $pnext) = @_;

    my $read_strand = is_reverse($flag) ? "R" : "F";
    my $mate_strand = is_mate_reverse($flag) ? "R" : "F";

    my $read_no = is_read1($flag) ? "1" : is_read2($flag) ? "2" : "?";
    my $mate_no;
    if ($read_no eq "1") {
        $mate_no = "2";
    } elsif ($read_no eq "2") {
        $mate_no = "1";
    } else {
        $mate_no = "?";
    }

    my $read_label = $read_strand . $read_no;
    my $mate_label = $mate_strand . $mate_no;

    if ($pos <= $pnext) {
        return $read_label . $mate_label;
    } else {
        return $mate_label . $read_label;
    }
}

sub median {
    my @x = sort { $a <=> $b } @_;
    return 0 unless @x;
    my $n = scalar @x;
    if ($n % 2) {
        return $x[int($n / 2)];
    } else {
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

sub log_msg {
    my ($verbose, $msg) = @_;
    print $msg . "\n" if $verbose;
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
  --config    Config file
  --bam       Coordinate-sorted BAM file
  --sample    Sample ID
  --out       Output discordant-read summary file

Config requirements based on conf/hcm_exondel.example.conf:
  SAMTOOLS
  REFSEQ_MANE_SELECT_GENE_TXT      required when SCAN_WHOLE_BAM=0
  HCM_CORE_GENE_LIST              optional; if set, only these genes are scanned when SCAN_WHOLE_BAM=0
  REF_FASTA_INDEX                  optional but recommended
  SCAN_WHOLE_BAM                  default: 0
  TARGET_REGION_FLANK             default: 0
  MIN_DISCORDANT_MAPQ             fallback: MIN_MAPQ, then 20
  MIN_DISCORDANT_INSERT_SIZE      default: 1000
  MIN_DISCORDANT_READS            default: 3
  DISCORDANT_CLUSTER_DISTANCE     default: 1000
  EXCLUDE_DUPLICATES              default: 1
  FILTER_DELETION_ORIENTATION    default: 1
  VERBOSE                         default: 1

Path rule:
  Relative config paths containing /, for example db/xxx, are resolved against
  the project root inferred from the config path.

Output:
  *.discordant_reads.tsv
  *.discordant_reads.raw_discordant_pairs.tsv
  *.discordant_reads.supporting_pairs.tsv
  *.discordant_reads.discarded_clusters.tsv

Notes:
  This script intentionally does not use MIN_DISCORDANT_CANDIDATE_SIZE
  or MAX_DISCORDANT_CANDIDATE_SIZE. Candidate-size filtering should be
  handled later during evidence merging or exon annotation.
USAGE
}



