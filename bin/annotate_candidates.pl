#!/usr/bin/env perl
# -*- coding: utf-8 -*-
#
# ============================================================
# Program:
#   annotate_candidates.pl
#
# Purpose:
#   Annotate high-confidence exon-level deletion candidates using
#   MANE RefSeq exon annotation.
#
# Main features:
#   1. Read final candidates from merged_candidates.tsv.
#   2. Read MANE exon annotation from REFSEQ_MANE_SELECT_EXON_TXT
#      configured in hcm_exondel.example.conf.
#   3. Annotate each candidate with affected gene, transcript and exon
#      information.
#   4. Write one compact review-friendly TSV file specified by --out.
#
#
# Important:
#
#   This script is designed for the new merge_evidence.pl output.
#
#   The input merged_candidates.tsv is assumed to contain final
#   high-confidence candidates only, for example candidates with:
#
#     Evidence_Level = High
#     Evidence_Count = 3
#     Evidence_Types = Depth,Split,Discordant
#
#   Therefore, this script no longer uses:
#
#     ANNOTATE_EVIDENCE_LEVEL
#
#   It no longer splits output into:
#
#     SAMPLE.Evidence_Level.High.tsv
#     SAMPLE.Evidence_Level.Other.tsv
#
#   All input rows are annotated and written directly to --out.
#
#
# Coordinate selection logic:
#
#   The new merged_candidates.tsv contains three coordinate systems:
#
#     1) Best_Start / Best_End
#
#        Recommended candidate boundary for downstream annotation and
#        reporting. This is selected by merge_evidence.pl using priority:
#
#          Split > Discordant > Depth
#
#        This is the preferred coordinate for exon annotation.
#
#
#     2) Merged_Start / Merged_End
#
#        Outer boundary of all overlapping evidence intervals.
#        This is usually the widest possible event range.
#        It is used only when Best_Start / Best_End are unavailable.
#
#
#     3) Core_Start / Core_End
#
#        Conservative overlap region supported by evidence types.
#        For three-evidence candidates, this is the overlap of
#        Depth + Split + Discordant.
#        It may be shorter than the real deletion.
#        It is used only when both Best and Merged coordinates are unavailable.
#
#
#   Coordinate priority:
#
#     Priority 1:
#       Chrom + Best_Start + Best_End
#
#     Priority 2:
#       Chrom + Merged_Start + Merged_End
#
#     Priority 3:
#       Chrom + Core_Start + Core_End
#
#   The selected coordinate source is written to output column:
#
#     Coordinate_Source
#
#   The output columns Start, End, Candidate_Region and Size_bp are based
#   on the selected coordinate system.
#
#
# Input:
#   --config  hcm_exondel config file
#   --input   merged_candidates.tsv
#   --out     annotated_candidates.tsv
#
# Required input columns:
#   First column: sample name
#   Chrom
#
# Required coordinate columns:
#   At least one of the following coordinate sets must exist:
#
#     Best_Start + Best_End
#     Merged_Start + Merged_End
#     Core_Start + Core_End
#
# Config parameters used:
#   REFSEQ_MANE_SELECT_EXON_TXT
#   MIN_EXON_OVERLAP_FRACTION
#
# Expected exon TXT format:
#   Gene    Transcript    Exon    Chrom    Start    End    Strand
#
# Coordinate system:
#   1-based closed interval
#
# Example:
#   perl bin/annotate_candidates.pl \
#     --config conf/hcm_exondel.example.conf \
#     --input  test/test_results/25B09089386/04.candidates/25B09089386.merged_candidates.tsv \
#     --out    test/test_results/25B09089386/06.report/25B09089386.annotated_candidates.tsv
# ============================================================

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path getcwd);

# ============================================================
# Parse command-line arguments
# ============================================================

my $config;
my $input;
my $out;
my $help = 0;

GetOptions(
    "config=s" => \$config,
    "input=s"  => \$input,
    "out=s"    => \$out,
    "help"     => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless $config && $input && $out;

$config = abs_path($config);
$input  = abs_path($input);

die "[ERROR] Config file not found: $config\n"
    unless defined $config && -s $config;

die "[ERROR] Input candidate file not found: $input\n"
    unless defined $input && -s $input;

# ============================================================
# Read config
# ============================================================

my %CONF = read_config($config);

my $project_root = guess_project_root($config);
my $config_dir   = dirname($config);

my $exon_txt = get_conf_required(\%CONF, "REFSEQ_MANE_SELECT_EXON_TXT");
$exon_txt = resolve_path($exon_txt, $project_root, $config_dir);

die "[ERROR] REFSEQ_MANE_SELECT_EXON_TXT file not found: $exon_txt\n"
    unless -s $exon_txt;

my $min_exon_overlap_fraction = get_conf_value(
    \%CONF,
    "MIN_EXON_OVERLAP_FRACTION",
    0.20
);

die "[ERROR] MIN_EXON_OVERLAP_FRACTION must be numeric: $min_exon_overlap_fraction\n"
    unless is_number($min_exon_overlap_fraction);

die "[ERROR] MIN_EXON_OVERLAP_FRACTION must be > 0 and <= 1: $min_exon_overlap_fraction\n"
    unless $min_exon_overlap_fraction > 0 && $min_exon_overlap_fraction <= 1;

# ============================================================
# Step 1. Read candidate file
# ============================================================

my ($header_ref, $rows_ref) = read_tsv($input);
my @header = @$header_ref;
my @rows   = @$rows_ref;
my %idx    = header_index(@header);

die "[ERROR] No candidate rows found in input: $input\n" unless @rows;

my $sample_col = $header[0];

die "[ERROR] The first column of input file is empty\n"
    unless defined $sample_col && $sample_col ne "";

my $file_sample = get_first_column_sample_name(
    \@rows,
    $sample_col,
    \@header,
    \%idx
);

$file_sample = sanitize_filename($file_sample);

die "[ERROR] Required column 'Chrom' not found in input: $input\n"
    unless exists $idx{Chrom};

my @coord_sets = get_supported_coordinate_sets(\%idx);

die "[ERROR] No supported coordinate columns found in input: $input\n"
  . "        Required one of: Best_Start/Best_End, Merged_Start/Merged_End, Core_Start/Core_End\n"
    unless @coord_sets;

# ============================================================
# Step 2. Prepare output file
# ============================================================

my $outdir = dirname($out);
$outdir = resolve_output_dir($outdir);
make_path($outdir) unless -d $outdir;

print "[INFO] Candidate annotation started\n";
print "[INFO] Config                    : $config\n";
print "[INFO] Project root              : $project_root\n";
print "[INFO] Input candidate file       : $input\n";
print "[INFO] Sample column             : $sample_col\n";
print "[INFO] Output file sample prefix : $file_sample\n";
print "[INFO] Exon TXT                  : $exon_txt\n";
print "[INFO] Min exon overlap fraction : $min_exon_overlap_fraction\n";
print "[INFO] Coordinate priority       : ", join(" > ", map { $_->{mode} } @coord_sets), "\n";
print "[INFO] Output file               : $out\n";

# ============================================================
# Step 3. Read MANE exon annotation
# ============================================================

my @exons = read_exon_txt($exon_txt);

die "[ERROR] No valid exon records found in REFSEQ_MANE_SELECT_EXON_TXT: $exon_txt\n"
    unless @exons;

my %exons_by_chr;

foreach my $exon (@exons) {
    foreach my $key (chrom_keys($exon->{chrom})) {
        push @{ $exons_by_chr{$key} }, $exon;
    }
}

foreach my $chr (keys %exons_by_chr) {
    @{ $exons_by_chr{$chr} } = sort {
        $a->{start} <=> $b->{start} ||
        $a->{end}   <=> $b->{end}
    } @{ $exons_by_chr{$chr} };
}

# ============================================================
# Step 4. Open output file
# ============================================================

my @compact_header = qw(
    Sample
    Chrom
    Start
    End
    Candidate_Region
    Size_bp
    Coordinate_Source
    Scan_Source
    Evidence_Level
    Evidence_Types
    Evidence_Count
    Depth_Support
    Split_Read_Support
    Discordant_Read_Support
    Annotated_Gene
    Annotated_Transcript
    Affected_Exons
    Overlap_Exon_Count
    Fully_Covered_Exons
    Partially_Overlapped_Exons
    Annotation_Status
    Exon_Overlap_Detail
);

open my $out_fh, ">", $out
    or die "[ERROR] Cannot write output file: $out\n";

print $out_fh join("\t", @compact_header) . "\n";

# ============================================================
# Step 5. Annotate all candidates
# ============================================================

my $total_count     = 0;
my $annotated_count = 0;
my $not_annotated   = 0;

foreach my $row (@rows) {
    $total_count++;

    my %r = row_hash($row, \@header, \%idx);

    my ($coord_source, $chr, $start, $end) = select_candidate_coordinate(
        \%r,
        \@coord_sets,
        $total_count
    );

    $r{__ANNOT_CHROM}  = $chr;
    $r{__ANNOT_START}  = $start;
    $r{__ANNOT_END}    = $end;
    $r{__COORD_SOURCE} = $coord_source;

    my $annotation = annotate_one_candidate(
        chrom                     => $chr,
        start                     => $start,
        end                       => $end,
        exons_by_chr              => \%exons_by_chr,
        min_exon_overlap_fraction => $min_exon_overlap_fraction,
    );

    if ($annotation->{annotation_status} eq "Annotated") {
        $annotated_count++;
    }
    else {
        $not_annotated++;
    }

    print_compact_output_row($out_fh, \%r, $annotation, $sample_col);
}

close $out_fh;

print "[INFO] Candidate annotation finished\n";
print "[INFO] Total candidates     : $total_count\n";
print "[INFO] Annotated candidates : $annotated_count\n";
print "[INFO] Not annotated        : $not_annotated\n";
print "[INFO] Output file          : $out\n";

exit 0;

# ============================================================
# Coordinate selection functions
# ============================================================

sub get_supported_coordinate_sets {
    my ($idx_ref) = @_;

    my @all = (
        {
            mode  => "Best",
            chrom => "Chrom",
            start => "Best_Start",
            end   => "Best_End",
        },
        {
            mode  => "Merged",
            chrom => "Chrom",
            start => "Merged_Start",
            end   => "Merged_End",
        },
        {
            mode  => "Core",
            chrom => "Chrom",
            start => "Core_Start",
            end   => "Core_End",
        },
    );

    my @supported;

    foreach my $set (@all) {
        next unless exists $idx_ref->{ $set->{chrom} };
        next unless exists $idx_ref->{ $set->{start} };
        next unless exists $idx_ref->{ $set->{end} };

        push @supported, $set;
    }

    return @supported;
}


sub select_candidate_coordinate {
    my ($row_hash_ref, $coord_sets_ref, $row_no) = @_;

    foreach my $set (@$coord_sets_ref) {
        my $chr   = get_value($row_hash_ref, $set->{chrom});
        my $start = get_value($row_hash_ref, $set->{start});
        my $end   = get_value($row_hash_ref, $set->{end});

        next if !defined $chr || $chr eq "" || $chr eq "NA";
        next unless is_integer($start);
        next unless is_integer($end);

        if ($start > $end) {
            ($start, $end) = ($end, $start);
        }

        return ($set->{mode}, $chr, $start, $end);
    }

    die "[ERROR] No valid coordinate found at candidate row $row_no\n";
}

# ============================================================
# Annotation functions
# ============================================================

sub annotate_one_candidate {
    my %args = @_;

    my $chr          = $args{chrom};
    my $start        = $args{start};
    my $end          = $args{end};
    my $exons_by_chr = $args{exons_by_chr};
    my $min_fraction = $args{min_exon_overlap_fraction};

    my @overlapped_exons;
    my @fully_covered;
    my @partial;
    my @details;

    my %genes;
    my %transcripts;

    if (!exists $exons_by_chr->{$chr}) {
        return empty_annotation("No exon annotation on chromosome $chr");
    }

    foreach my $exon (@{ $exons_by_chr->{$chr} }) {
        last if $exon->{start} > $end;
        next if $exon->{end} < $start;

        my $overlap_len = overlap_length(
            $start,
            $end,
            $exon->{start},
            $exon->{end}
        );

        next if $overlap_len <= 0;

        my $exon_len = $exon->{end} - $exon->{start} + 1;
        my $frac     = $overlap_len / $exon_len;

        next if $frac < $min_fraction;

        my $exon_label = join(
            "|",
            $exon->{gene},
            $exon->{transcript},
            "EX" . $exon->{exon}
        );

        push @overlapped_exons, $exon_label;

        $genes{ $exon->{gene} }++;
        $transcripts{ $exon->{transcript} }++;

        my $detail = join(
            ":",
            $exon_label,
            $exon->{chrom},
            $exon->{start},
            $exon->{end},
            "overlap=" . $overlap_len,
            "fraction=" . sprintf("%.3f", $frac)
        );

        push @details, $detail;

        if ($start <= $exon->{start} && $end >= $exon->{end}) {
            push @fully_covered, $exon_label;
        }
        else {
            push @partial, $exon_label;
        }
    }

    if (!@overlapped_exons) {
        return empty_annotation("No exon overlapped");
    }

    my $gene       = choose_most_frequent_key(%genes)       || "NA";
    my $transcript = choose_most_frequent_key(%transcripts) || "NA";

    my @uniq_overlapped = unique(@overlapped_exons);

    return {
        gene                       => $gene,
        transcript                 => $transcript,
        affected_exons             => join(",", @uniq_overlapped),
        overlap_exon_count         => scalar(@uniq_overlapped),
        fully_covered_exons        => @fully_covered ? join(",", unique(@fully_covered)) : "NA",
        partially_overlapped_exons => @partial ? join(",", unique(@partial)) : "NA",
        exon_overlap_detail        => join(";", @details),
        annotation_status          => "Annotated",
    };
}


sub empty_annotation {
    my ($reason) = @_;

    return {
        gene                       => "NA",
        transcript                 => "NA",
        affected_exons             => "NA",
        overlap_exon_count         => 0,
        fully_covered_exons        => "NA",
        partially_overlapped_exons => "NA",
        exon_overlap_detail        => $reason,
        annotation_status          => "Not_annotated",
    };
}


sub overlap_length {
    my ($s1, $e1, $s2, $e2) = @_;

    my $s = $s1 > $s2 ? $s1 : $s2;
    my $e = $e1 < $e2 ? $e1 : $e2;

    return 0 if $e < $s;
    return $e - $s + 1;
}

# ============================================================
# Read MANE RefSeq exon TXT
# ============================================================

sub read_exon_txt {
    my ($file) = @_;

    my @exons;

    open my $fh, "<", $file
        or die "[ERROR] Cannot open exon TXT: $file\n";

    my $header = <$fh>;
    die "[ERROR] Empty exon TXT file: $file\n" unless defined $header;

    chomp $header;
    $header =~ s/\r$//;

    my @header = split /\t/, $header, -1;
    my %idx    = header_index(@header);

    foreach my $required (qw/Gene Transcript Exon Chrom Start End Strand/) {
        die "[ERROR] Required column '$required' not found in exon TXT: $file\n"
            unless exists $idx{$required};
    }

    my $line_no = 1;

    while (my $line = <$fh>) {
        $line_no++;

        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        my @f = split /\t/, $line, -1;

        my $gene       = get_field(\@f, \%idx, "Gene");
        my $transcript = get_field(\@f, \%idx, "Transcript");
        my $exon       = get_field(\@f, \%idx, "Exon");
        my $chrom      = get_field(\@f, \%idx, "Chrom");
        my $start      = get_field(\@f, \%idx, "Start");
        my $end        = get_field(\@f, \%idx, "End");
        my $strand     = get_field(\@f, \%idx, "Strand");

        die "[ERROR] Invalid exon coordinate at $file line $line_no: $line\n"
            unless is_integer($start) && is_integer($end) && $start <= $end;

        $gene       = "NA" if !defined $gene       || $gene eq "";
        $transcript = "NA" if !defined $transcript || $transcript eq "";
        $exon       = "NA" if !defined $exon       || $exon eq "";
        $chrom      = "NA" if !defined $chrom      || $chrom eq "";
        $strand     = "NA" if !defined $strand     || $strand eq "";

        $exon =~ s/^EX//i;

        push @exons, {
            gene       => $gene,
            transcript => $transcript,
            exon       => $exon,
            chrom      => $chrom,
            start      => $start,
            end        => $end,
            strand     => $strand,
        };
    }

    close $fh;

    return @exons;
}

# ============================================================
# Compact output helper
# ============================================================

sub print_compact_output_row {
    my ($fh, $row_hash_ref, $annotation, $sample_col) = @_;

    my $sample_name = get_value($row_hash_ref, $sample_col);
    $sample_name = "NA" if !defined $sample_name || $sample_name eq "";

    my $chrom = get_value($row_hash_ref, "__ANNOT_CHROM");
    my $start = get_value($row_hash_ref, "__ANNOT_START");
    my $end   = get_value($row_hash_ref, "__ANNOT_END");

    my $coordinate_source = get_value($row_hash_ref, "__COORD_SOURCE");

    my $candidate_region = "NA";
    if ($chrom ne "NA" && is_integer($start) && is_integer($end)) {
        $candidate_region = "$chrom:$start-$end";
    }

    my $size_bp = "NA";
    if (is_integer($start) && is_integer($end)) {
        my $s = $start;
        my $e = $end;
        ($s, $e) = ($e, $s) if $s > $e;
        $size_bp = $e - $s + 1;
    }

    my $scan_source = first_existing_value(
        $row_hash_ref,
        qw/Scan_Source Candidate_Source Source/
    );

    if ($coordinate_source ne "NA") {
        my $best_evidence = first_existing_value(
            $row_hash_ref,
            qw/Best_Evidence/
        );

        if ($best_evidence ne "NA" && $coordinate_source eq "Best") {
            $scan_source = "$coordinate_source;Best_Evidence=$best_evidence";
        }
        elsif ($scan_source eq "NA") {
            $scan_source = $coordinate_source;
        }
        else {
            $scan_source = "$coordinate_source;$scan_source";
        }
    }

    my $evidence_level = get_value($row_hash_ref, "Evidence_Level");

    my $evidence_types = first_existing_value(
        $row_hash_ref,
        qw/Evidence_Types Evidence_Type Evidence_Sources Support_Types Support_Type Evidence/
    );

    my $evidence_count = first_existing_value(
        $row_hash_ref,
        qw/Evidence_Count Support_Count Num_Evidence Evidence_Type_Count/
    );

    my $depth_support = first_existing_value(
        $row_hash_ref,
        qw/Depth_Support Has_Depth Depth_Evidence Depth_Candidate Depth/
    );

    my $split_read_support = first_existing_value(
        $row_hash_ref,
        qw/Split_Read_Support Has_Split_Read Split_Evidence Split_Read_Count Split_Reads Split_Read/
    );

    my $discordant_read_support = first_existing_value(
        $row_hash_ref,
        qw/Discordant_Read_Support Has_Discordant_Read Discordant_Evidence Discordant_Read_Count Discordant_Reads Discordant_Pair_Support Discordant_Pair/
    );

    my @out = (
        $sample_name,
        $chrom,
        $start,
        $end,
        $candidate_region,
        $size_bp,
        $coordinate_source,
        $scan_source,
        $evidence_level,
        $evidence_types,
        $evidence_count,
        $depth_support,
        $split_read_support,
        $discordant_read_support,
        $annotation->{gene},
        $annotation->{transcript},
        $annotation->{affected_exons},
        $annotation->{overlap_exon_count},
        $annotation->{fully_covered_exons},
        $annotation->{partially_overlapped_exons},
        $annotation->{annotation_status},
        $annotation->{exon_overlap_detail},
    );

    print $fh join("\t", @out) . "\n";
}

# ============================================================
# TSV helpers
# ============================================================

sub read_tsv {
    my ($file) = @_;

    open my $fh, "<", $file
        or die "[ERROR] Cannot open file: $file\n";

    my $header = <$fh>;
    die "[ERROR] Empty input file: $file\n" unless defined $header;

    chomp $header;
    $header =~ s/\r$//;

    my @header = split /\t/, $header, -1;
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
        $idx{ $header[$i] } = $i;
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


sub get_field {
    my ($fields_ref, $idx_ref, $key) = @_;

    my $i = $idx_ref->{$key};

    return defined $fields_ref->[$i] ? $fields_ref->[$i] : "";
}


sub get_value {
    my ($row_hash_ref, $key) = @_;

    if (
        exists $row_hash_ref->{$key}
        && defined $row_hash_ref->{$key}
        && $row_hash_ref->{$key} ne ""
    ) {
        return $row_hash_ref->{$key};
    }

    return "NA";
}


sub first_existing_value {
    my ($row_hash_ref, @keys) = @_;

    foreach my $key (@keys) {
        if (
            exists $row_hash_ref->{$key}
            && defined $row_hash_ref->{$key}
            && $row_hash_ref->{$key} ne ""
            && $row_hash_ref->{$key} ne "NA"
        ) {
            return $row_hash_ref->{$key};
        }
    }

    return "NA";
}


sub get_first_column_sample_name {
    my ($rows_ref, $sample_col, $header_ref, $idx_ref) = @_;

    die "[ERROR] No candidate rows found in input file\n"
        unless @$rows_ref;

    my %seen_sample;
    my $first_sample = "";

    foreach my $row (@$rows_ref) {
        my %r = row_hash($row, $header_ref, $idx_ref);

        my $s = get_value(\%r, $sample_col);

        next if $s eq "NA" || $s eq "";

        $seen_sample{$s}++;

        if ($first_sample eq "") {
            $first_sample = $s;
        }
    }

    die "[ERROR] Cannot determine sample name from the first column: $sample_col\n"
        if $first_sample eq "";

    if (scalar(keys %seen_sample) > 1) {
        print STDERR "[WARN] Multiple sample names found in first column '$sample_col'. ";
        print STDERR "Output file prefix will use the first sample: $first_sample\n";
    }

    return $first_sample;
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

        $line =~ s/\s+#.*$//;

        if ($line =~ /^\s*([^=\s]+)\s*=\s*(.*?)\s*$/) {
            my $key = $1;
            my $val = $2;

            $val =~ s/^\s+//;
            $val =~ s/\s+$//;
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

# ============================================================
# Path helpers
# ============================================================

sub guess_project_root {
    my ($config_file) = @_;

    my $dir = dirname($config_file);

    if ($dir =~ /\/conf$/) {
        $dir =~ s/\/conf$//;
        return $dir;
    }

    return dirname($config_file);
}


sub resolve_path {
    my ($path, $project_root, $config_dir) = @_;

    die "[ERROR] Empty path\n" unless defined $path && $path ne "";

    if ($path =~ /^\//) {
        return abs_path($path) || $path;
    }

    my $p1 = "$project_root/$path";
    return abs_path($p1) if -e $p1;

    my $p2 = "$config_dir/$path";
    return abs_path($p2) if -e $p2;

    my $p3 = getcwd() . "/$path";
    return abs_path($p3) if -e $p3;

    return $p1;
}


sub resolve_output_dir {
    my ($dir) = @_;

    die "[ERROR] Empty output directory\n" unless defined $dir && $dir ne "";

    if ($dir =~ /^\//) {
        return $dir;
    }

    return getcwd() . "/$dir";
}

# ============================================================
# Chromosome helpers
# ============================================================

sub chrom_keys {
    my ($chrom) = @_;

    return ("NA") if !defined $chrom || $chrom eq "";

    my @keys;
    my %seen;

    push @keys, $chrom;
    $seen{$chrom} = 1;

    if ($chrom =~ /^chr(.+)$/i) {
        my $no_chr = $1;
        push @keys, $no_chr unless $seen{$no_chr}++;
    }
    else {
        my $with_chr = "chr$chrom";
        push @keys, $with_chr unless $seen{$with_chr}++;
    }

    return @keys;
}

# ============================================================
# Misc helpers
# ============================================================

sub sanitize_filename {
    my ($x) = @_;

    $x = "NA" if !defined $x || $x eq "";

    $x =~ s/^\s+//;
    $x =~ s/\s+$//;
    $x =~ s/[\/\\:\*\?\"\<\>\|\s]+/_/g;

    return $x;
}


sub choose_most_frequent_key {
    my (%hash) = @_;

    return "" unless %hash;

    my @sorted = sort {
        $hash{$b} <=> $hash{$a} ||
        $a cmp $b
    } keys %hash;

    return $sorted[0];
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


sub is_integer {
    my ($x) = @_;

    return defined $x && $x =~ /^\d+$/;
}


sub is_number {
    my ($x) = @_;

    return defined $x && $x =~ /^-?(?:\d+\.?\d*|\.\d+)$/;
}


sub usage {
    return <<"USAGE";

Usage:
  perl bin/annotate_candidates.pl \\
    --config conf/hcm_exondel.example.conf \\
    --input  test_results/SAMPLE/04.candidates/SAMPLE.merged_candidates.tsv \\
    --out    test_results/SAMPLE/06.report/SAMPLE.annotated_candidates.tsv

Purpose:
  Annotate all candidates in merged_candidates.tsv.

Note:
  The script does not use ANNOTATE_EVIDENCE_LEVEL.
  It does not split High/Other outputs.
  All input rows are annotated and written directly to --out.

Coordinate priority:
  1. Chrom + Best_Start + Best_End
  2. Chrom + Merged_Start + Merged_End
  3. Chrom + Core_Start + Core_End

Config parameters used:
  REFSEQ_MANE_SELECT_EXON_TXT
  MIN_EXON_OVERLAP_FRACTION

Compact output columns:
  Sample
  Chrom
  Start
  End
  Candidate_Region
  Size_bp
  Coordinate_Source
  Scan_Source
  Evidence_Level
  Evidence_Types
  Evidence_Count
  Depth_Support
  Split_Read_Support
  Discordant_Read_Support
  Annotated_Gene
  Annotated_Transcript
  Affected_Exons
  Overlap_Exon_Count
  Fully_Covered_Exons
  Partially_Overlapped_Exons
  Annotation_Status
  Exon_Overlap_Detail

Expected exon TXT format:
  Gene    Transcript    Exon    Chrom    Start    End    Strand

USAGE
}


