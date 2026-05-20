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
#   1. depth reduction
#   2. discordant read pairs
#   3. split / soft-clipped reads
#
# Main logic:
#   - Merge records by chromosome and genomic proximity.
#   - Evidence records are merged if:
#       record.start <= current_cluster.end + EVIDENCE_MERGE_DISTANCE
#   - Candidate must be supported by at least two evidence types.
#
# Config:
#   EVIDENCE_MERGE_DISTANCE=1000
#   MIN_EVIDENCE_COUNT=2
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

$config = safe_abs_path($config);
die "[ERROR] Config file not found or empty: $config\n" unless defined $config && -s $config;

$depth_file      = safe_abs_path($depth_file);
$discordant_file = safe_abs_path($discordant_file);
$split_file      = safe_abs_path($split_file);

my %CONF = read_config($config);

my $merge_distance = get_conf_value(\%CONF, "EVIDENCE_MERGE_DISTANCE", 1000);
my $min_evidence_count = get_conf_value(\%CONF, "MIN_EVIDENCE_COUNT", 2);

die "[ERROR] EVIDENCE_MERGE_DISTANCE must be a non-negative integer. Current value: $merge_distance\n"
    unless defined $merge_distance && $merge_distance =~ /^\d+$/;

die "[ERROR] MIN_EVIDENCE_COUNT must be 2 or 3. Current value: $min_evidence_count\n"
    unless defined $min_evidence_count
        && $min_evidence_count =~ /^\d+$/
        && $min_evidence_count >= 2
        && $min_evidence_count <= 3;

my $outdir = dirname($out);
make_path($outdir) unless -d $outdir;

print "[INFO] Merge evidence started\n";
print "[INFO] Sample             : $sample\n";
print "[INFO] Config             : $config\n";
print "[INFO] Depth file         : " . defined_or_na($depth_file) . "\n";
print "[INFO] Discordant file    : " . defined_or_na($discordant_file) . "\n";
print "[INFO] Split file         : " . defined_or_na($split_file) . "\n";
print "[INFO] Merge distance     : $merge_distance\n";
print "[INFO] Min evidence count : $min_evidence_count\n";

# ============================================================
# Step 1. Read evidence files
# ============================================================

my @depth_records = ();
my @discordant_records = ();
my @split_records = ();

if (file_has_content($depth_file)) {
    @depth_records = read_depth_evidence($depth_file);
} else {
    warn "[WARN] Depth evidence file missing or empty: " . defined_or_na($depth_file) . "\n";
}

if (file_has_content($discordant_file)) {
    @discordant_records = read_discordant_evidence($discordant_file);
} else {
    warn "[WARN] Discordant evidence file missing or empty: " . defined_or_na($discordant_file) . "\n";
}

if (file_has_content($split_file)) {
    @split_records = read_split_evidence($split_file);
} else {
    warn "[WARN] Split evidence file missing or empty: " . defined_or_na($split_file) . "\n";
}

my @all_records = (
    @depth_records,
    @discordant_records,
    @split_records,
);

# ============================================================
# Step 2. Open output and print header
# ============================================================

open my $out_fh, ">", $out or die "[ERROR] Cannot write output file: $out\n";

print $out_fh join("\t", qw(
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
    Min_Mean_Depth_Ratio
    Max_Mean_Depth_Ratio
    Discordant_Reads
    Median_Insert_Size
    Split_Reads
    Left_Softclip_Reads
    Right_Softclip_Reads
    SA_Tag_Reads
    Depth_Record_Count
    Discordant_Record_Count
    Split_Record_Count
    Source_Record_Count
    Source_Records
    Source_Targets
    All_Genes
    All_Transcripts
    All_Regions
    Candidate_Status
    Comment
)) . "\n";

if (!@all_records) {
    close $out_fh;
    print "[INFO] No evidence records found. Empty report generated: $out\n";
    exit 0;
}

# ============================================================
# Step 3. Cluster evidence records
# ============================================================

my @clusters = cluster_evidence_records(
    records        => \@all_records,
    merge_distance => $merge_distance,
);

# ============================================================
# Step 4. Build final candidates
# ============================================================

my $pass_count = 0;

foreach my $cluster (@clusters) {
    my $candidate = build_final_candidate($cluster);

    next if $candidate->{evidence_count} < $min_evidence_count;

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
        $candidate->{has_depth} ? "Yes" : "No",
        $candidate->{has_discordant} ? "Yes" : "No",
        $candidate->{has_split} ? "Yes" : "No",
        $candidate->{depth_status},
        $candidate->{region_mean_depth},
        $candidate->{candidate_mean_depth},
        $candidate->{mean_depth_ratio},
        $candidate->{min_mean_depth_ratio},
        $candidate->{max_mean_depth_ratio},
        $candidate->{discordant_reads},
        $candidate->{median_insert_size},
        $candidate->{split_reads},
        $candidate->{left_softclip_reads},
        $candidate->{right_softclip_reads},
        $candidate->{sa_tag_reads},
        $candidate->{depth_record_count},
        $candidate->{discordant_record_count},
        $candidate->{split_record_count},
        $candidate->{source_record_count},
        $candidate->{source_records},
        $candidate->{source_targets},
        $candidate->{all_genes},
        $candidate->{all_transcripts},
        $candidate->{all_regions},
        $candidate->{candidate_status},
        $candidate->{comment},
    ) . "\n";
}

close $out_fh;

print "[INFO] Merge evidence finished\n";
print "[INFO] Depth records      : " . scalar(@depth_records) . "\n";
print "[INFO] Discordant records : " . scalar(@discordant_records) . "\n";
print "[INFO] Split records      : " . scalar(@split_records) . "\n";
print "[INFO] Evidence clusters  : " . scalar(@clusters) . "\n";
print "[INFO] Passed candidates  : $pass_count\n";
print "[INFO] Final report       : $out\n";

exit 0;

# ============================================================
# Read depth evidence
# ============================================================

sub read_depth_evidence {
    my ($file) = @_;

    my @records;

    my ($header_ref, $rows_ref) = read_tsv($file);
    my @header = @$header_ref;
    my @rows   = @$rows_ref;

    return @records unless @header;

    my %idx = header_index(@header);

    # Region is not required because current depth_candidates.tsv does not have this column.
    # Candidate_Length is also optional and can be calculated from Start/End.
    foreach my $required (qw(SampleID Gene Transcript Chrom Start End)) {
        die "[ERROR] Required column '$required' not found in depth file: $file\n"
            unless exists $idx{$required};
    }

    foreach my $row (@rows) {
        my %r = row_hash($row, \@header, \%idx);

        next unless valid_coordinate($r{Chrom}, $r{Start}, $r{End});

        my $candidate_status = exists $idx{Candidate_Status}
            ? ($r{Candidate_Status} || "NA")
            : "Depth_candidate";

        next if $candidate_status ne "Depth_candidate";

        my $gene       = clean_value($r{Gene});
        my $transcript = clean_value($r{Transcript});

        my $region = "NA";
        if (exists $idx{Region}) {
            $region = clean_value($r{Region});
        }
        elsif (exists $idx{Exon}) {
            $region = clean_value($r{Exon});
        }
        elsif (exists $idx{Target_Region}) {
            $region = clean_value($r{Target_Region});
        }
        elsif (exists $idx{Target_Name}) {
            my ($g2, $t2, $r2) = parse_target_name($r{Target_Name});
            $region = $r2;
        }

        my $candidate_length = exists $idx{Candidate_Length}
            ? ($r{Candidate_Length} || ($r{End} - $r{Start} + 1))
            : ($r{End} - $r{Start} + 1);

        push @records, {
            evidence_type        => "Depth",
            sample               => clean_value($r{SampleID}),
            gene                 => $gene,
            transcript           => $transcript,
            region               => $region,
            chr                  => clean_value($r{Chrom}),
            start                => int($r{Start}),
            end                  => int($r{End}),
            length               => $candidate_length,
            depth_status         => exists $idx{Depth_Status} ? clean_value($r{Depth_Status}) : "NA",
            region_mean_depth    => exists $idx{Region_Mean_Depth} ? clean_value($r{Region_Mean_Depth}) : "NA",
            candidate_mean_depth => exists $idx{Candidate_Mean_Depth} ? clean_value($r{Candidate_Mean_Depth}) : "NA",
            mean_depth_ratio     => exists $idx{Mean_Depth_Ratio} ? clean_value($r{Mean_Depth_Ratio}) : "NA",
            target_name          => join("|", $gene, $transcript, $region),
            source_record        => join(":", "Depth", $gene, $transcript, $region, $r{Chrom}, $r{Start}, $r{End}),
        };
    }

    return @records;
}


# ============================================================
# Read discordant evidence
# ============================================================

sub read_discordant_evidence {
    my ($file) = @_;

    my @records;

    my ($header_ref, $rows_ref) = read_tsv($file);
    my @header = @$header_ref;
    my @rows   = @$rows_ref;

    return @records unless @header;

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

        my $target_name = clean_value($r{Target_Name} || "NA");
        my ($gene, $transcript, $region) = parse_target_name($target_name);

        push @records, {
            evidence_type      => "Discordant",
            sample             => clean_value($r{SampleID}),
            gene               => $gene,
            transcript         => $transcript,
            region             => $region,
            chr                => clean_value($r{Chrom}),
            start              => int($r{Start}),
            end                => int($r{End}),
            length             => $r{Candidate_Length} || ($r{End} - $r{Start} + 1),
            discordant_reads   => clean_value($r{Discordant_Reads} || 0),
            median_insert_size => clean_value($r{Median_Insert_Size} || "NA"),
            target_name        => $target_name,
            source_record      => join(":", "Discordant", $gene, $transcript, $region, $r{Chrom}, $r{Start}, $r{End}),
        };
    }

    return @records;
}

# ============================================================
# Read split evidence
# ============================================================

sub read_split_evidence {
    my ($file) = @_;

    my @records;

    my ($header_ref, $rows_ref) = read_tsv($file);
    my @header = @$header_ref;
    my @rows   = @$rows_ref;

    return @records unless @header;

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

        my $target_name = clean_value($r{Target_Name} || "NA");
        my ($gene, $transcript, $region) = parse_target_name($target_name);

        push @records, {
            evidence_type        => "Split",
            sample               => clean_value($r{SampleID}),
            gene                 => $gene,
            transcript           => $transcript,
            region               => $region,
            chr                  => clean_value($r{Chrom}),
            start                => int($r{Start}),
            end                  => int($r{End}),
            length               => $r{Candidate_Length} || ($r{End} - $r{Start} + 1),
            split_reads          => clean_value($r{Split_Reads} || 0),
            left_softclip_reads  => clean_value($r{Left_Softclip_Reads} || 0),
            right_softclip_reads => clean_value($r{Right_Softclip_Reads} || 0),
            sa_tag_reads         => clean_value($r{SA_Tag_Reads} || 0),
            target_name          => $target_name,
            source_record        => join(":", "Split", $gene, $transcript, $region, $r{Chrom}, $r{Start}, $r{End}),
        };
    }

    return @records;
}

# ============================================================
# Cluster evidence records
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
        my $current_end;

        foreach my $r (@sorted) {
            if (!@current) {
                @current = ($r);
                $current_end = $r->{end};
                next;
            }

            if ($r->{start} <= $current_end + $merge_distance) {
                push @current, $r;
                $current_end = max($current_end, $r->{end});
            } else {
                push @clusters, [@current];

                @current = ($r);
                $current_end = $r->{end};
            }
        }

        push @clusters, [@current] if @current;
    }

    return @clusters;
}

# ============================================================
# Build final candidate
# ============================================================

sub build_final_candidate {
    my ($cluster_ref) = @_;

    my @records = @$cluster_ref;

    my $chr = $records[0]->{chr};

    my @starts = map { $_->{start} } @records;
    my @ends   = map { $_->{end} } @records;

    my $start = min(@starts);
    my $end   = max(@ends);

    my %evidence_types;
    my %genes;
    my %transcripts;
    my %regions;
    my %targets;

    my $has_depth      = 0;
    my $has_discordant = 0;
    my $has_split      = 0;

    my @depth_status;
    my @region_mean_depth;
    my @candidate_mean_depth;
    my @mean_depth_ratio;

    my $discordant_reads = 0;
    my @median_insert_size;

    my $split_reads          = 0;
    my $left_softclip_reads  = 0;
    my $right_softclip_reads = 0;
    my $sa_tag_reads         = 0;

    my $depth_record_count      = 0;
    my $discordant_record_count = 0;
    my $split_record_count      = 0;

    my @source_records;

    foreach my $r (@records) {
        $evidence_types{ $r->{evidence_type} } = 1;

        $genes{ $r->{gene} }++
            if is_known_value($r->{gene});

        $transcripts{ $r->{transcript} }++
            if is_known_value($r->{transcript});

        $regions{ $r->{region} }++
            if is_known_value($r->{region});

        $targets{ $r->{target_name} }++
            if is_known_value($r->{target_name});

        push @source_records, $r->{source_record};

        if ($r->{evidence_type} eq "Depth") {
            $has_depth = 1;
            $depth_record_count++;

            push @depth_status, $r->{depth_status}
                if is_known_value($r->{depth_status});

            push @region_mean_depth, $r->{region_mean_depth}
                if is_numeric($r->{region_mean_depth});

            push @candidate_mean_depth, $r->{candidate_mean_depth}
                if is_numeric($r->{candidate_mean_depth});

            push @mean_depth_ratio, $r->{mean_depth_ratio}
                if is_numeric($r->{mean_depth_ratio});
        }
        elsif ($r->{evidence_type} eq "Discordant") {
            $has_discordant = 1;
            $discordant_record_count++;

            $discordant_reads += numeric_or_zero($r->{discordant_reads});

            push @median_insert_size, $r->{median_insert_size}
                if is_numeric($r->{median_insert_size});
        }
        elsif ($r->{evidence_type} eq "Split") {
            $has_split = 1;
            $split_record_count++;

            $split_reads          += numeric_or_zero($r->{split_reads});
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

    my @unique_source_records = unique(@source_records);

    my $comment = build_comment(
        has_depth      => $has_depth,
        has_discordant => $has_discordant,
        has_split      => $has_split,
    );

    return {
        gene             => $gene,
        transcript       => $transcript,
        region           => $region,
        chr              => $chr,
        start            => $start,
        end              => $end,
        candidate_length => $end - $start + 1,

        evidence_count => $evidence_count,
        evidence_types => $evidence_types,
        evidence_level => $evidence_level,

        has_depth      => $has_depth,
        has_discordant => $has_discordant,
        has_split      => $has_split,

        depth_status => @depth_status ? join(",", unique(@depth_status)) : "NA",

        region_mean_depth => @region_mean_depth
            ? sprintf("%.4f", mean(@region_mean_depth))
            : "NA",

        candidate_mean_depth => @candidate_mean_depth
            ? sprintf("%.4f", mean(@candidate_mean_depth))
            : "NA",

        mean_depth_ratio => @mean_depth_ratio
            ? sprintf("%.4f", mean(@mean_depth_ratio))
            : "NA",

        min_mean_depth_ratio => @mean_depth_ratio
            ? sprintf("%.4f", min(@mean_depth_ratio))
            : "NA",

        max_mean_depth_ratio => @mean_depth_ratio
            ? sprintf("%.4f", max(@mean_depth_ratio))
            : "NA",

        discordant_reads => $discordant_reads,

        median_insert_size => @median_insert_size
            ? sprintf("%.2f", mean(@median_insert_size))
            : "NA",

        split_reads          => $split_reads,
        left_softclip_reads  => $left_softclip_reads,
        right_softclip_reads => $right_softclip_reads,
        sa_tag_reads         => $sa_tag_reads,

        depth_record_count      => $depth_record_count,
        discordant_record_count => $discordant_record_count,
        split_record_count      => $split_record_count,
        source_record_count     => scalar(@unique_source_records),

        source_records => join(",", @unique_source_records),
        source_targets => %targets ? join(",", sort keys %targets) : "NA",

        all_genes       => %genes ? join(",", sort keys %genes) : "NA",
        all_transcripts => %transcripts ? join(",", sort keys %transcripts) : "NA",
        all_regions     => %regions ? join(",", sort keys %regions) : "NA",

        candidate_status => "Pass",
        comment          => $comment,
    };
}

# ============================================================
# Classification and comments
# ============================================================

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
# TSV helpers
# ============================================================

sub read_tsv {
    my ($file) = @_;

    open my $fh, "<", $file or die "[ERROR] Cannot open file: $file\n";

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

# ============================================================
# Config helpers
# ============================================================

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

    if (
        exists $conf_ref->{$key}
        && defined $conf_ref->{$key}
        && $conf_ref->{$key} ne ""
    ) {
        return $conf_ref->{$key};
    }

    return $default;
}

# ============================================================
# General helpers
# ============================================================

sub parse_target_name {
    my ($target_name) = @_;

    return ("NA", "NA", "NA")
        unless defined $target_name && $target_name ne "";

    my @targets = split /,/, $target_name;
    my $name = "NA";

    foreach my $t (@targets) {
        $t = clean_value($t);

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

    return (
        clean_value($gene),
        clean_value($transcript),
        clean_value($region),
    );
}

sub valid_coordinate {
    my ($chr, $start, $end) = @_;

    return 0 unless defined $chr && $chr ne "" && $chr ne "NA";
    return 0 unless defined $start && $start =~ /^\d+$/;
    return 0 unless defined $end   && $end   =~ /^\d+$/;
    return 0 unless $end >= $start;

    return 1;
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
    return $x if is_numeric($x);

    return 0;
}

sub is_numeric {
    my ($x) = @_;

    return 0 unless defined $x;
    return 0 if $x eq "NA";
    return $x =~ /^-?\d+(\.\d+)?$/ ? 1 : 0;
}

sub mean {
    my @x = @_;

    my $sum = 0;
    my $n   = 0;

    foreach my $v (@x) {
        next unless is_numeric($v);

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

        $v = clean_value($v);

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

sub clean_value {
    my ($x) = @_;

    return "NA" unless defined $x;

    $x =~ s/^\s+//;
    $x =~ s/\s+$//;

    return "NA" if $x eq "";

    return $x;
}

sub is_known_value {
    my ($x) = @_;

    return 0 unless defined $x;

    $x = clean_value($x);

    return 0 if $x eq "";
    return 0 if $x eq "NA";
    return 0 if $x eq ".";
    return 0 if $x eq "whole_bam";

    return 1;
}

sub file_has_content {
    my ($file) = @_;

    return 0 unless defined $file;
    return 0 unless -e $file;
    return 0 unless -s $file;

    return 1;
}

sub safe_abs_path {
    my ($path) = @_;

    return undef unless defined $path && $path ne "";

    my $abs = abs_path($path);

    return defined $abs ? $abs : $path;
}

sub defined_or_na {
    my ($x) = @_;

    return defined $x ? $x : "NA";
}

sub usage {
    return <<"USAGE";

Usage:
  perl bin/merge_evidence.pl \\
    --config conf/hcm_exondel.conf \\
    --sample SAMPLE001 \\
    --depth results/SAMPLE001/01.depth/SAMPLE001.depth_candidates.tsv \\
    --discordant results/SAMPLE001/03.discordant_reads/SAMPLE001.discordant_reads.tsv \\
    --split results/SAMPLE001/02.split_reads/SAMPLE001.split_reads.tsv \\
    --out results/SAMPLE001/04.candidates/SAMPLE001.merged_candidates.tsv

Required arguments:
  --config       Config file
  --sample       Sample ID
  --depth        Depth candidate file
  --discordant   Discordant read candidate file
  --split        Split-read candidate file
  --out          Final merged candidate report

Config parameters used by this script:
  EVIDENCE_MERGE_DISTANCE
  MIN_EVIDENCE_COUNT

Recommended config:
  EVIDENCE_MERGE_DISTANCE=1000
  MIN_EVIDENCE_COUNT=2

Output:
  merged_candidates.tsv

USAGE
}

