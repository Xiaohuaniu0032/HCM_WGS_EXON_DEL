#!/usr/bin/env perl
# -*- coding: utf-8 -*-

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path);

# ============================================================
# merge_evidence.pl
#
# Merge evidence from:
#   1. window-based depth analysis
#   2. discordant read-pair analysis
#   3. split-read / soft-clipped read analysis
#
# Design principle:
#   Evidence records are merged by genomic overlap or proximity,
#   not by exact Start/End equality.
#
# Required command-line arguments:
#   --config
#   --sample
#   --depth
#   --discordant
#   --split
#   --out
#
# Config-controlled parameters:
#   EVIDENCE_MERGE_DISTANCE
#   MIN_EVIDENCE_COUNT
#   REQUIRE_DEPTH_EVIDENCE
# ============================================================

my $config;
my $sample;
my $depth_file;
my $discordant_file;
my $split_file;
my $out;
my $help = 0;

GetOptions(
    "config=s"     => \$config,
    "sample=s"     => \$sample,
    "depth=s"      => \$depth_file,
    "discordant=s" => \$discordant_file,
    "split=s"      => \$split_file,
    "out=s"        => \$out,
    "help"         => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless $config && $sample && $depth_file && $discordant_file && $split_file && $out;

$config = abs_path($config);

die "[ERROR] Config file not found: $config\n" unless -s $config;

$depth_file      = abs_path($depth_file);
$discordant_file = abs_path($discordant_file);
$split_file      = abs_path($split_file);

die "[ERROR] Depth evidence file not found: $depth_file\n" unless -s $depth_file;
die "[ERROR] Discordant evidence file not found: $discordant_file\n" unless -s $discordant_file;
die "[ERROR] Split evidence file not found: $split_file\n" unless -s $split_file;

my %CONF = read_config($config);

my $merge_distance = get_conf_value(\%CONF, "EVIDENCE_MERGE_DISTANCE", 1000);
my $min_evidence_count = get_conf_value(\%CONF, "MIN_EVIDENCE_COUNT", 2);
my $require_depth_evidence = normalize_bool(
    get_conf_value(\%CONF, "REQUIRE_DEPTH_EVIDENCE", 0)
);

my $outdir = dirname($out);
make_path($outdir) unless -d $outdir;

print "[INFO] Merge evidence started\n";
print "[INFO] Sample                 : $sample\n";
print "[INFO] Depth file             : $depth_file\n";
print "[INFO] Discordant file        : $discordant_file\n";
print "[INFO] Split file             : $split_file\n";
print "[INFO] Merge distance         : $merge_distance\n";
print "[INFO] Min evidence count     : $min_evidence_count\n";
print "[INFO] Require depth evidence : $require_depth_evidence\n";

# ============================================================
# Step 1. Read three evidence files
# ============================================================

my @depth_records      = read_depth_evidence($depth_file);
my @discordant_records = read_discordant_evidence($discordant_file);
my @split_records      = read_split_evidence($split_file);

my @all_records = (
    @depth_records,
    @discordant_records,
    @split_records,
);

# If no evidence, still output header.
open my $out_fh, ">", $out
    or die "[ERROR] Cannot write output file: $out\n";

print $out_fh join("\t",
    qw(
        SampleID
        Gene
        Transcript
        Region
        Chrom
        Start
        End
        Candidate_Length
        Evidence_Count
        Evidence_Types
        Evidence_Level
        Depth_Evidence
        Discordant_Evidence
        Split_Evidence
        Depth_Status
        Region_Mean_Depth
        Candidate_Mean_Depth
        Mean_Depth_Ratio
        Discordant_Reads
        Median_Insert_Size
        Split_Reads
        Left_Softclip_Reads
        Right_Softclip_Reads
        SA_Tag_Reads
        Source_Records
        Candidate_Status
        Comment
    )
) . "\n";

if (!@all_records) {
    close $out_fh;
    print "[INFO] No evidence records found. Empty report generated: $out\n";
    exit 0;
}

# ============================================================
# Step 2. Cluster evidence records by genomic coordinate
# ============================================================

my @clusters = cluster_evidence_records(
    records        => \@all_records,
    merge_distance => $merge_distance,
);

# ============================================================
# Step 3. Build final candidates
# ============================================================

my $pass_count = 0;

foreach my $cluster (@clusters) {
    my $candidate = build_final_candidate($cluster);

    next if $candidate->{evidence_count} < $min_evidence_count;

    if ($require_depth_evidence && !$candidate->{has_depth}) {
        next;
    }

    $pass_count++;

    print $out_fh join("\t",
        $sample,
        $candidate->{gene},
        $candidate->{transcript},
        $candidate->{region},
        $candidate->{chr},
        $candidate->{start},
        $candidate->{end},
        $candidate->{candidate_length},
        $candidate->{evidence_count},
        $candidate->{evidence_types},
        $candidate->{evidence_level},
        $candidate->{has_depth}      ? "Yes" : "No",
        $candidate->{has_discordant} ? "Yes" : "No",
        $candidate->{has_split}      ? "Yes" : "No",
        $candidate->{depth_status},
        $candidate->{region_mean_depth},
        $candidate->{candidate_mean_depth},
        $candidate->{mean_depth_ratio},
        $candidate->{discordant_reads},
        $candidate->{median_insert_size},
        $candidate->{split_reads},
        $candidate->{left_softclip_reads},
        $candidate->{right_softclip_reads},
        $candidate->{sa_tag_reads},
        $candidate->{source_records},
        $candidate->{candidate_status},
        $candidate->{comment},
    ) . "\n";
}

close $out_fh;

print "[INFO] Merge evidence finished\n";
print "[INFO] Depth records       : " . scalar(@depth_records) . "\n";
print "[INFO] Discordant records  : " . scalar(@discordant_records) . "\n";
print "[INFO] Split records       : " . scalar(@split_records) . "\n";
print "[INFO] Evidence clusters   : " . scalar(@clusters) . "\n";
print "[INFO] Passed candidates   : $pass_count\n";
print "[INFO] Final report        : $out\n";

exit 0;


# ============================================================
# Read evidence files
# ============================================================

sub read_depth_evidence {
    my ($file) = @_;

    my @records;
    my ($header_ref, $rows_ref) = read_tsv($file);

    my @header = @$header_ref;
    my @rows   = @$rows_ref;

    my %idx = header_index(@header);

    foreach my $required (qw(SampleID Gene Transcript Region Chrom Start End Candidate_Length)) {
        die "[ERROR] Required column '$required' not found in depth file: $file\n"
            unless exists $idx{$required};
    }

    foreach my $row (@rows) {
        my %r = row_hash($row, \@header, \%idx);

        next unless valid_coordinate($r{Chrom}, $r{Start}, $r{End});

        my $depth_status = $r{Depth_Status} || "NA";
        my $candidate_status = $r{Candidate_Status} || "NA";

        next if $candidate_status ne "Depth_candidate";

        push @records, {
            evidence_type          => "Depth",
            sample                 => $r{SampleID},
            gene                   => $r{Gene} || "NA",
            transcript             => $r{Transcript} || "NA",
            region                 => $r{Region} || "NA",
            chr                    => $r{Chrom},
            start                  => $r{Start},
            end                    => $r{End},
            length                 => $r{Candidate_Length} || ($r{End} - $r{Start} + 1),
            depth_status           => $depth_status,
            region_mean_depth      => $r{Region_Mean_Depth} || "NA",
            candidate_mean_depth   => $r{Candidate_Mean_Depth} || "NA",
            mean_depth_ratio       => $r{Mean_Depth_Ratio} || "NA",
            source_record          => join(":", "Depth", $r{Chrom}, $r{Start}, $r{End}),
        };
    }

    return @records;
}


sub read_discordant_evidence {
    my ($file) = @_;

    my @records;
    my ($header_ref, $rows_ref) = read_tsv($file);

    my @header = @$header_ref;
    my @rows   = @$rows_ref;

    my %idx = header_index(@header);

    foreach my $required (qw(SampleID Chrom Start End Candidate_Length Discordant_Reads)) {
        die "[ERROR] Required column '$required' not found in discordant file: $file\n"
            unless exists $idx{$required};
    }

    foreach my $row (@rows) {
        my %r = row_hash($row, \@header, \%idx);

        next unless valid_coordinate($r{Chrom}, $r{Start}, $r{End});

        my $candidate_status = $r{Candidate_Status} || "NA";
        next if $candidate_status ne "Discordant_candidate";

        my $target_name = $r{Target_Name} || "NA";
        my ($gene, $transcript, $region) = parse_target_name($target_name);

        push @records, {
            evidence_type        => "Discordant",
            sample               => $r{SampleID},
            gene                 => $gene,
            transcript           => $transcript,
            region               => $region,
            chr                  => $r{Chrom},
            start                => $r{Start},
            end                  => $r{End},
            length               => $r{Candidate_Length} || ($r{End} - $r{Start} + 1),
            discordant_reads     => $r{Discordant_Reads} || 0,
            median_insert_size   => $r{Median_Insert_Size} || "NA",
            target_name          => $target_name,
            source_record        => join(":", "Discordant", $r{Chrom}, $r{Start}, $r{End}),
        };
    }

    return @records;
}


sub read_split_evidence {
    my ($file) = @_;

    my @records;
    my ($header_ref, $rows_ref) = read_tsv($file);

    my @header = @$header_ref;
    my @rows   = @$rows_ref;

    my %idx = header_index(@header);

    foreach my $required (qw(SampleID Chrom Start End Candidate_Length Split_Reads)) {
        die "[ERROR] Required column '$required' not found in split file: $file\n"
            unless exists $idx{$required};
    }

    foreach my $row (@rows) {
        my %r = row_hash($row, \@header, \%idx);

        next unless valid_coordinate($r{Chrom}, $r{Start}, $r{End});

        my $candidate_status = $r{Candidate_Status} || "NA";
        next if $candidate_status ne "Split_candidate";

        my $target_name = $r{Target_Name} || "NA";
        my ($gene, $transcript, $region) = parse_target_name($target_name);

        push @records, {
            evidence_type          => "Split",
            sample                 => $r{SampleID},
            gene                   => $gene,
            transcript             => $transcript,
            region                 => $region,
            chr                    => $r{Chrom},
            start                  => $r{Start},
            end                    => $r{End},
            length                 => $r{Candidate_Length} || ($r{End} - $r{Start} + 1),
            split_reads            => $r{Split_Reads} || 0,
            left_softclip_reads    => $r{Left_Softclip_Reads} || 0,
            right_softclip_reads   => $r{Right_Softclip_Reads} || 0,
            sa_tag_reads           => $r{SA_Tag_Reads} || 0,
            target_name            => $target_name,
            source_record          => join(":", "Split", $r{Chrom}, $r{Start}, $r{End}),
        };
    }

    return @records;
}


# ============================================================
# Cluster and build candidates
# ============================================================

sub cluster_evidence_records {
    my %args = @_;

    my $records_ref    = $args{records};
    my $merge_distance = $args{merge_distance};

    my @records = @$records_ref;

    my %by_chr;

    foreach my $r (@records) {
        push @{ $by_chr{ $r->{chr} } }, $r;
    }

    my @clusters;

    foreach my $chr (sort keys %by_chr) {
        my @sorted = sort {
            $a->{start} <=> $b->{start}
            ||
            $a->{end} <=> $b->{end}
        } @{ $by_chr{$chr} };

        my @current;
        my $current_start;
        my $current_end;

        foreach my $r (@sorted) {
            if (!@current) {
                @current = ($r);
                $current_start = $r->{start};
                $current_end   = $r->{end};
                next;
            }

            if ($r->{start} <= $current_end + $merge_distance) {
                push @current, $r;
                $current_start = $r->{start} if $r->{start} < $current_start;
                $current_end   = $r->{end}   if $r->{end}   > $current_end;
            }
            else {
                push @clusters, [@current];
                @current = ($r);
                $current_start = $r->{start};
                $current_end   = $r->{end};
            }
        }

        if (@current) {
            push @clusters, [@current];
        }
    }

    return @clusters;
}


sub build_final_candidate {
    my ($cluster_ref) = @_;

    my @records = @$cluster_ref;

    my $chr = $records[0]->{chr};

    my @starts = map { $_->{start} } @records;
    my @ends   = map { $_->{end}   } @records;

    my $start = min(@starts);
    my $end   = max(@ends);

    my %evidence_types;
    my %genes;
    my %transcripts;
    my %regions;

    my $has_depth      = 0;
    my $has_discordant = 0;
    my $has_split      = 0;

    my @depth_status;
    my @region_mean_depth;
    my @candidate_mean_depth;
    my @mean_depth_ratio;

    my $discordant_reads = 0;
    my @median_insert_size;

    my $split_reads = 0;
    my $left_softclip_reads = 0;
    my $right_softclip_reads = 0;
    my $sa_tag_reads = 0;

    my @source_records;

    foreach my $r (@records) {
        $evidence_types{ $r->{evidence_type} } = 1;

        $genes{ $r->{gene} }++ if defined $r->{gene} && $r->{gene} ne "" && $r->{gene} ne "NA";
        $transcripts{ $r->{transcript} }++ if defined $r->{transcript} && $r->{transcript} ne "" && $r->{transcript} ne "NA";
        $regions{ $r->{region} }++ if defined $r->{region} && $r->{region} ne "" && $r->{region} ne "NA";

        push @source_records, $r->{source_record};

        if ($r->{evidence_type} eq "Depth") {
            $has_depth = 1;

            push @depth_status, $r->{depth_status} if defined $r->{depth_status};
            push @region_mean_depth, $r->{region_mean_depth}
                if defined $r->{region_mean_depth} && $r->{region_mean_depth} ne "NA";
            push @candidate_mean_depth, $r->{candidate_mean_depth}
                if defined $r->{candidate_mean_depth} && $r->{candidate_mean_depth} ne "NA";
            push @mean_depth_ratio, $r->{mean_depth_ratio}
                if defined $r->{mean_depth_ratio} && $r->{mean_depth_ratio} ne "NA";
        }
        elsif ($r->{evidence_type} eq "Discordant") {
            $has_discordant = 1;

            $discordant_reads += numeric_or_zero($r->{discordant_reads});

            push @median_insert_size, $r->{median_insert_size}
                if defined $r->{median_insert_size} && $r->{median_insert_size} ne "NA";
        }
        elsif ($r->{evidence_type} eq "Split") {
            $has_split = 1;

            $split_reads += numeric_or_zero($r->{split_reads});
            $left_softclip_reads  += numeric_or_zero($r->{left_softclip_reads});
            $right_softclip_reads += numeric_or_zero($r->{right_softclip_reads});
            $sa_tag_reads         += numeric_or_zero($r->{sa_tag_reads});
        }
    }

    my $evidence_count = scalar(keys %evidence_types);
    my $evidence_types = join(",", sort keys %evidence_types);

    my $evidence_level = classify_evidence_level($evidence_count);

    my $gene       = choose_most_frequent_key(%genes);
    my $transcript = choose_most_frequent_key(%transcripts);
    my $region     = choose_most_frequent_key(%regions);

    $gene       = "NA" unless $gene;
    $transcript = "NA" unless $transcript;
    $region     = "NA" unless $region;

    my $candidate_status = $evidence_count >= 2 ? "Pass" : "Low_confidence";

    my $comment = build_comment(
        has_depth      => $has_depth,
        has_discordant => $has_discordant,
        has_split      => $has_split,
        evidence_count => $evidence_count,
    );

    return {
        gene                   => $gene,
        transcript             => $transcript,
        region                 => $region,
        chr                    => $chr,
        start                  => $start,
        end                    => $end,
        candidate_length       => $end - $start + 1,
        evidence_count         => $evidence_count,
        evidence_types         => $evidence_types,
        evidence_level         => $evidence_level,
        has_depth              => $has_depth,
        has_discordant         => $has_discordant,
        has_split              => $has_split,
        depth_status           => @depth_status ? join(",", unique(@depth_status)) : "NA",
        region_mean_depth      => @region_mean_depth ? sprintf("%.4f", mean(@region_mean_depth)) : "NA",
        candidate_mean_depth   => @candidate_mean_depth ? sprintf("%.4f", mean(@candidate_mean_depth)) : "NA",
        mean_depth_ratio       => @mean_depth_ratio ? sprintf("%.4f", mean(@mean_depth_ratio)) : "NA",
        discordant_reads       => $discordant_reads,
        median_insert_size     => @median_insert_size ? sprintf("%.2f", mean(@median_insert_size)) : "NA",
        split_reads            => $split_reads,
        left_softclip_reads    => $left_softclip_reads,
        right_softclip_reads   => $right_softclip_reads,
        sa_tag_reads           => $sa_tag_reads,
        source_records         => join(",", unique(@source_records)),
        candidate_status       => $candidate_status,
        comment                => $comment,
    };
}


sub classify_evidence_level {
    my ($count) = @_;

    return "High"   if $count >= 3;
    return "Medium" if $count == 2;
    return "Low"    if $count == 1;

    return "None";
}


sub build_comment {
    my %args = @_;

    my $has_depth      = $args{has_depth};
    my $has_discordant = $args{has_discordant};
    my $has_split      = $args{has_split};
    my $evidence_count = $args{evidence_count};

    if ($has_depth && $has_discordant && $has_split) {
        return "Supported by depth reduction, discordant read pairs, and split reads";
    }

    if ($has_depth && $has_discordant) {
        return "Supported by depth reduction and discordant read pairs";
    }

    if ($has_depth && $has_split) {
        return "Supported by depth reduction and split reads";
    }

    if ($has_discordant && $has_split) {
        return "Supported by discordant read pairs and split reads";
    }

    if ($has_depth) {
        return "Supported only by depth reduction";
    }

    if ($has_discordant) {
        return "Supported only by discordant read pairs";
    }

    if ($has_split) {
        return "Supported only by split reads";
    }

    return "No valid evidence";
}


# ============================================================
# Generic TSV helpers
# ============================================================

sub read_tsv {
    my ($file) = @_;

    open my $fh, "<", $file
        or die "[ERROR] Cannot open file: $file\n";

    my $header = <$fh>;

    if (!defined $header) {
        close $fh;
        return ([], []);
    }

    chomp $header;
    $header =~ s/\r$//;

    my @header = split /\t/, $header;

    my @rows;

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;

        my @f = split /\t/, $line, -1;
        push @rows, \@f;
    }

    close $fh;

    return (\@header, \@rows);
}


sub header_index {
    my (@header) = @_;

    my %idx;

    for (my $i = 0; $i < @header; $i++) {
        $idx{$header[$i]} = $i;
    }

    return %idx;
}


sub row_hash {
    my ($row_ref, $header_ref, $idx_ref) = @_;

    my @row    = @$row_ref;
    my @header = @$header_ref;

    my %r;

    foreach my $h (@header) {
        my $i = $idx_ref->{$h};
        $r{$h} = defined $row[$i] ? $row[$i] : "NA";
    }

    return %r;
}


sub valid_coordinate {
    my ($chr, $start, $end) = @_;

    return 0 unless defined $chr && $chr ne "" && $chr ne "NA";
    return 0 unless defined $start && $start =~ /^\d+$/;
    return 0 unless defined $end && $end =~ /^\d+$/;
    return 0 unless $end >= $start;

    return 1;
}


# ============================================================
# Config helpers
# ============================================================

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


# ============================================================
# Misc helpers
# ============================================================

sub parse_target_name {
    my ($target_name) = @_;

    return ("NA", "NA", "NA") unless defined $target_name && $target_name ne "";

    # If multiple target names are merged, use the first non-NA target.
    my @targets = split /,/, $target_name;
    my $name = "NA";

    foreach my $t (@targets) {
        if ($t && $t ne "NA" && $t ne "whole_bam") {
            $name = $t;
            last;
        }
    }

    $name = $targets[0] if $name eq "NA" && @targets;

    my ($gene, $transcript, $region);

    if ($name =~ /\|/) {
        my @x = split /\|/, $name;
        $gene       = $x[0] || "NA";
        $transcript = $x[1] || "NA";
        $region     = $x[2] || "NA";
    }
    elsif ($name =~ /:/) {
        my @x = split /:/, $name;
        $gene       = $x[0] || "NA";
        $transcript = $x[1] || "NA";
        $region     = $x[2] || "NA";
    }
    else {
        $gene       = $name || "NA";
        $transcript = "NA";
        $region     = "NA";
    }

    return ($gene, $transcript, $region);
}


sub choose_most_frequent_key {
    my (%hash) = @_;

    return "" unless %hash;

    my @sorted = sort {
        $hash{$b} <=> $hash{$a}
        ||
        $a cmp $b
    } keys %hash;

    return $sorted[0];
}


sub numeric_or_zero {
    my ($x) = @_;

    return 0 unless defined $x;
    return 0 if $x eq "NA";
    return $x if $x =~ /^-?\d+(\.\d+)?$/;

    return 0;
}


sub mean {
    my @x = @_;

    return 0 unless @x;

    my $sum = 0;
    my $n = 0;

    foreach my $v (@x) {
        next unless defined $v;
        next if $v eq "NA";
        next unless $v =~ /^-?\d+(\.\d+)?$/;

        $sum += $v;
        $n++;
    }

    return 0 unless $n > 0;

    return $sum / $n;
}


sub unique {
    my @x = @_;

    my %seen;
    my @u;

    foreach my $v (@x) {
        next unless defined $v;
        next if $v eq "";

        next if $seen{$v}++;
        push @u, $v;
    }

    return @u;
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


sub usage {
    return <<"USAGE";

Usage:

  perl bin/merge_evidence.pl \\
    --config conf/hcm_exondel.conf \\
    --sample SAMPLE001 \\
    --depth results/SAMPLE001/01.depth/SAMPLE001.depth_candidates.tsv \\
    --discordant results/SAMPLE001/02.discordant_reads/SAMPLE001.discordant_reads.tsv \\
    --split results/SAMPLE001/03.split_reads/SAMPLE001.split_reads.tsv \\
    --out results/SAMPLE001/05.report/SAMPLE001.final_report.tsv

Required arguments:

  --config       Config file
  --sample       Sample ID
  --depth        Depth candidate file
  --discordant   Discordant read candidate file
  --split        Split-read candidate file
  --out          Final merged report

Config parameters used by this script:

  EVIDENCE_MERGE_DISTANCE
  MIN_EVIDENCE_COUNT
  REQUIRE_DEPTH_EVIDENCE

Recommended config:

  EVIDENCE_MERGE_DISTANCE=1000
  MIN_EVIDENCE_COUNT=2
  REQUIRE_DEPTH_EVIDENCE=0

Output:

  final_report.tsv

USAGE
}

