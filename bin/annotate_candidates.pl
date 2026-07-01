#!/usr/bin/env perl
# -*- coding: utf-8 -*-

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path getcwd);

# ============================================================
# Program:
#   annotate_candidates.pl
#
# Purpose:
#   Annotate merged exon-level deletion candidates using MANE RefSeq
#   exon annotation.
#
# Current workflow design:
#   merge_evidence.pl is split-read-centered:
#     1. Split-read cluster is used as the candidate event backbone.
#     2. Depth evidence and discordant-read evidence are used as
#        supporting validation.
#
#   Therefore, this script annotates all rows from merged_candidates.tsv,
#   including:
#     High
#     Moderate
#     Low
#
#   It does not filter by Evidence_Level.
#   It does not split output into High/Other files.
#
# Important:
#   If the input merged_candidates.tsv contains only a header and no
#   candidate rows, this script writes an output file with only the
#   annotated header and exits normally.
#
# Coordinate selection logic:
#   Priority 1:
#     Chrom + Best_Start + Best_End
#
#   Priority 2:
#     Chrom + Split_Start + Split_End
#
#   Priority 3:
#     Chrom + Merged_Start + Merged_End
#
#   Priority 4:
#     Chrom + Core_Start + Core_End
#
#   The selected coordinate source is written to:
#     Coordinate_Source
#
# Input:
#   --config hcm_exondel config file
#   --input  merged_candidates.tsv
#   --out    annotated_candidates.tsv
#
# Config parameters used:
#   REFSEQ_MANE_SELECT_EXON_TXT
#   MIN_EXON_OVERLAP_FRACTION
#
# Expected exon TXT format:
#   Gene Transcript Exon Chrom Start End Strand
#
# Coordinate system:
#   1-based closed interval
#
# Example:
#   perl bin/annotate_candidates.pl \
#     --config conf/hcm_exondel.example.conf \
#     --input test_results/SAMPLE/04.candidates/SAMPLE.merged_candidates.tsv \
#     --out test_results/SAMPLE/05.report/SAMPLE.annotated_candidates.tsv
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

die usage() unless defined $config && defined $input && defined $out;

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
# Read candidate file
# ============================================================

my ($header_ref, $rows_ref) = read_tsv($input);

my @header = @$header_ref;
my @rows   = @$rows_ref;
my %idx    = header_index(@header);

die "[ERROR] Input candidate file has no header columns: $input\n"
    unless @header;

my $sample_col = $header[0];

die "[ERROR] The first column of input file is empty\n"
    unless defined $sample_col && $sample_col ne "";

die "[ERROR] Required column 'Chrom' not found in input: $input\n"
    unless exists $idx{Chrom};

my @coord_sets = get_supported_coordinate_sets(\%idx);

die "[ERROR] No supported coordinate columns found in input: $input\n"
  . "        Required one of:\n"
  . "        Best_Start/Best_End\n"
  . "        Split_Start/Split_End\n"
  . "        Merged_Start/Merged_End\n"
  . "        Core_Start/Core_End\n"
    unless @coord_sets;

# ============================================================
# Prepare output
# ============================================================

my @annotated_header = annotated_output_header();

my $outdir = dirname($out);
$outdir = resolve_output_dir($outdir);
make_path($outdir) unless -d $outdir;

my $file_sample = "NA";

if (@rows) {
    $file_sample = get_first_column_sample_name(
        \@rows,
        $sample_col,
        \@header,
        \%idx
    );
    $file_sample = sanitize_filename($file_sample);
}

print "[INFO] Candidate annotation started\n";
print "[INFO] Config : $config\n";
print "[INFO] Project root : $project_root\n";
print "[INFO] Input candidate file : $input\n";
print "[INFO] Candidate rows : ", scalar(@rows), "\n";
print "[INFO] Sample column : $sample_col\n";
print "[INFO] Output file sample prefix : $file_sample\n";
print "[INFO] Exon TXT : $exon_txt\n";
print "[INFO] Min exon overlap fraction : $min_exon_overlap_fraction\n";
print "[INFO] Coordinate priority : ", join(" > ", map { $_->{mode} } @coord_sets), "\n";
print "[INFO] Output file : $out\n";

# If merged_candidates.tsv has only header, write header-only report.
if (!@rows) {
    open my $empty_fh, ">", $out
        or die "[ERROR] Cannot write output file: $out\n";

    print $empty_fh join("\t", @annotated_header), "\n";
    close $empty_fh;

    print "[INFO] No candidate rows found. Header-only annotation file generated.\n";
    print "[INFO] Output file : $out\n";
    exit 0;
}

# ============================================================
# Read MANE exon annotation
# ============================================================

my @exons = read_exon_txt($exon_txt);

die "[ERROR] No valid exon records found in REFSEQ_MANE_SELECT_EXON_TXT: $exon_txt\n"
    unless @exons;

my %exons_by_chr;

for my $exon (@exons) {
    for my $key (chrom_keys($exon->{chrom})) {
        push @{ $exons_by_chr{$key} }, $exon;
    }
}

for my $chr (keys %exons_by_chr) {
    @{ $exons_by_chr{$chr} } = sort {
           $a->{start} <=> $b->{start}
        || $a->{end}   <=> $b->{end}
    } @{ $exons_by_chr{$chr} };
}

# ============================================================
# Annotate candidates
# ============================================================

open my $out_fh, ">", $out
    or die "[ERROR] Cannot write output file: $out\n";

print $out_fh join("\t", @annotated_header), "\n";

my $total_count     = 0;
my $annotated_count = 0;
my $not_annotated   = 0;

for my $row (@rows) {
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

    print_annotated_output_row(
        $out_fh,
        \%r,
        $annotation,
        $sample_col
    );
}

close $out_fh;

print "[INFO] Candidate annotation finished\n";
print "[INFO] Total candidates : $total_count\n";
print "[INFO] Annotated candidates : $annotated_count\n";
print "[INFO] Not annotated : $not_annotated\n";
print "[INFO] Output file : $out\n";

exit 0;

# ============================================================
# Output header
# ============================================================

sub annotated_output_header {
    return qw(
        SampleID
        Gene
        Chrom
        Start
        End
        Candidate_Region
        Size_bp
        Coordinate_Source

        Cluster_ID
        Split_Start
        Split_End
        Split_Size
        Split_Read_Support
        Split_Read_Count
        Split_Record_Count

        Depth_Support
        Depth_Covered_Bases
        Depth_Coverage_Fraction
        Depth_Record_Count
        Depth_Range

        Discordant_Read_Support
        Discordant_Overlap_Bases
        Discordant_Overlap_Fraction
        Discordant_Record_Count
        Discordant_Range
        Discordant_Cluster_IDs
        Discordant_Reads
        Median_Insert_Size

        Evidence_Level
        Evidence_Types
        Evidence_Count
        Best_Evidence
        Candidate_Status
        Source_Records
        Comment

        Annotated_Gene
        Annotated_Transcript
        Affected_Exons
        Overlap_Exon_Count
        Fully_Covered_Exons
        Partially_Overlapped_Exons
        Annotation_Status
        Exon_Overlap_Detail
    );
}

# ============================================================
# Coordinate selection
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
            mode  => "Split",
            chrom => "Chrom",
            start => "Split_Start",
            end   => "Split_End",
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

    for my $set (@all) {
        next unless exists $idx_ref->{ $set->{chrom} };
        next unless exists $idx_ref->{ $set->{start} };
        next unless exists $idx_ref->{ $set->{end} };

        push @supported, $set;
    }

    return @supported;
}

sub select_candidate_coordinate {
    my ($row_hash_ref, $coord_sets_ref, $row_no) = @_;

    for my $set (@$coord_sets_ref) {
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
# Annotation
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

    for my $exon (@{ $exons_by_chr->{$chr} }) {
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
        gene                         => $gene,
        transcript                   => $transcript,
        affected_exons               => join(",", @uniq_overlapped),
        overlap_exon_count           => scalar(@uniq_overlapped),
        fully_covered_exons          => @fully_covered ? join(",", unique(@fully_covered)) : "NA",
        partially_overlapped_exons   => @partial       ? join(",", unique(@partial))       : "NA",
        exon_overlap_detail          => join(";", @details),
        annotation_status            => "Annotated",
    };
}

sub empty_annotation {
    my ($reason) = @_;

    return {
        gene                         => "NA",
        transcript                   => "NA",
        affected_exons               => "NA",
        overlap_exon_count           => 0,
        fully_covered_exons          => "NA",
        partially_overlapped_exons   => "NA",
        exon_overlap_detail          => $reason,
        annotation_status            => "Not_annotated",
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
# Read exon TXT
# ============================================================

sub read_exon_txt {
    my ($file) = @_;

    my @exons;

    open my $fh, "<", $file
        or die "[ERROR] Cannot open exon TXT: $file\n";

    my $header = <$fh>;

    die "[ERROR] Empty exon TXT file: $file\n"
        unless defined $header;

    chomp $header;
    $header =~ s/\r$//;

    my @header = split /\t/, $header, -1;
    my %idx = header_index(@header);

    for my $required (qw/Gene Transcript Exon Chrom Start End Strand/) {
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
# Output row
# ============================================================

sub print_annotated_output_row {
    my ($fh, $row_hash_ref, $annotation, $sample_col) = @_;

    my $sample_name = first_existing_value(
        $row_hash_ref,
        "SampleID",
        "Sample",
        $sample_col
    );

    my $input_gene = first_existing_value(
        $row_hash_ref,
        "Gene",
        "Scan_Genes",
        "Genes",
        "Target_Gene",
        "Target_Name"
    );

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

    my $evidence_types = first_existing_value(
        $row_hash_ref,
        qw/Evidence_Types Evidence_Type Evidence_Sources Support_Types Support_Type Evidence/
    );

    my $cluster_id = first_existing_value(
        $row_hash_ref,
        qw/Cluster_ID ID Candidate_ID/
    );

    my $split_start = first_existing_value(
        $row_hash_ref,
        qw/Split_Start Start/
    );

    my $split_end = first_existing_value(
        $row_hash_ref,
        qw/Split_End End/
    );

    my $split_size = first_existing_value(
        $row_hash_ref,
        qw/Split_Size Best_Size Size Size_bp/
    );

    my $split_read_count = first_existing_value(
        $row_hash_ref,
        qw/Split_Reads Split_Read_Count Support_Reads/
    );

    my $split_record_count = first_existing_value(
        $row_hash_ref,
        qw/Split_Records Split_Record_Count Support_Records/
    );

    my $split_read_support = infer_split_support(
        row_hash_ref      => $row_hash_ref,
        evidence_types    => $evidence_types,
        split_read_count  => $split_read_count,
        cluster_id        => $cluster_id,
    );

    my $depth_support = infer_support(
        row_hash_ref      => $row_hash_ref,
        explicit_keys     => [qw/Depth_Support Has_Depth Depth_Evidence Depth_Candidate Depth/],
        evidence_types    => $evidence_types,
        evidence_name     => "Depth",
        numeric_keys      => [qw/Depth_Covered_Bases Depth_Record_Count/],
    );

    my $discordant_read_support = infer_support(
        row_hash_ref      => $row_hash_ref,
        explicit_keys     => [qw/Discordant_Read_Support Has_Discordant_Read Discordant_Evidence Discordant_Pair_Support Discordant_Pair/],
        evidence_types    => $evidence_types,
        evidence_name     => "Discordant",
        numeric_keys      => [qw/Discordant_Reads Discordant_Record_Count Discordant_Overlap_Bases/],
    );

    my @out = (
        $sample_name,
        $input_gene,
        $chrom,
        $start,
        $end,
        $candidate_region,
        $size_bp,
        $coordinate_source,

        $cluster_id,
        $split_start,
        $split_end,
        $split_size,
        $split_read_support,
        $split_read_count,
        $split_record_count,

        $depth_support,
        first_existing_value($row_hash_ref, qw/Depth_Covered_Bases/),
        first_existing_value($row_hash_ref, qw/Depth_Coverage_Fraction/),
        first_existing_value($row_hash_ref, qw/Depth_Record_Count/),
        first_existing_value($row_hash_ref, qw/Depth_Range/),

        $discordant_read_support,
        first_existing_value($row_hash_ref, qw/Discordant_Overlap_Bases/),
        first_existing_value($row_hash_ref, qw/Discordant_Overlap_Fraction/),
        first_existing_value($row_hash_ref, qw/Discordant_Record_Count/),
        first_existing_value($row_hash_ref, qw/Discordant_Range/),
        first_existing_value($row_hash_ref, qw/Discordant_Cluster_IDs/),
        first_existing_value($row_hash_ref, qw/Discordant_Reads/),
        first_existing_value($row_hash_ref, qw/Median_Insert_Size/),

        get_value($row_hash_ref, "Evidence_Level"),
        $evidence_types,
        first_existing_value($row_hash_ref, qw/Evidence_Count Support_Count Num_Evidence Evidence_Type_Count/),
        first_existing_value($row_hash_ref, qw/Best_Evidence/),
        first_existing_value($row_hash_ref, qw/Candidate_Status Status/),
        first_existing_value($row_hash_ref, qw/Source_Records/),
        first_existing_value($row_hash_ref, qw/Comment/),

        $annotation->{gene},
        $annotation->{transcript},
        $annotation->{affected_exons},
        $annotation->{overlap_exon_count},
        $annotation->{fully_covered_exons},
        $annotation->{partially_overlapped_exons},
        $annotation->{annotation_status},
        $annotation->{exon_overlap_detail},
    );

    @out = map { clean_tsv_value($_) } @out;

    print $fh join("\t", @out), "\n";
}

sub infer_split_support {
    my %args = @_;

    my $row_hash_ref     = $args{row_hash_ref};
    my $evidence_types   = $args{evidence_types};
    my $split_read_count = $args{split_read_count};
    my $cluster_id       = $args{cluster_id};

    my $explicit = first_existing_value(
        $row_hash_ref,
        qw/Split_Read_Support Has_Split_Read Split_Evidence/
    );

    return normalize_support_value($explicit) if $explicit ne "NA";

    return "Yes" if defined $evidence_types && $evidence_types =~ /(?:^|,)Split(?:,|$)/i;

    return "Yes" if is_numeric_positive($split_read_count);

    return "Yes" if defined $cluster_id && $cluster_id ne "" && $cluster_id ne "NA";

    return "No";
}

sub infer_support {
    my %args = @_;

    my $row_hash_ref   = $args{row_hash_ref};
    my $explicit_keys  = $args{explicit_keys};
    my $evidence_types = $args{evidence_types};
    my $evidence_name  = $args{evidence_name};
    my $numeric_keys   = $args{numeric_keys};

    my $explicit = first_existing_value($row_hash_ref, @$explicit_keys);

    return normalize_support_value($explicit) if $explicit ne "NA";

    if (defined $evidence_types && $evidence_types ne "NA") {
        my @types = split /,/, $evidence_types;

        for my $t (@types) {
            $t =~ s/^\s+|\s+$//g;
            return "Yes" if lc($t) eq lc($evidence_name);
        }
    }

    for my $key (@$numeric_keys) {
        my $v = first_existing_value($row_hash_ref, $key);
        return "Yes" if is_numeric_positive($v);
    }

    return "No";
}

sub normalize_support_value {
    my ($v) = @_;

    return "NA" unless defined $v;
    $v =~ s/^\s+|\s+$//g;

    return "NA" if $v eq "" || $v eq "NA";

    return "Yes" if $v =~ /^(1|yes|true|pass|support|supported)$/i;
    return "No"  if $v =~ /^(0|no|false|fail|none|unsupported)$/i;

    return "Yes" if is_numeric_positive($v);

    return $v;
}

# ============================================================
# TSV helpers
# ============================================================

sub read_tsv {
    my ($file) = @_;

    open my $fh, "<", $file
        or die "[ERROR] Cannot open file: $file\n";

    my $header = <$fh>;

    die "[ERROR] Empty input file: $file\n"
        unless defined $header;

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
        die "[ERROR] Empty column name found in input header at column " . ($i + 1) . "\n"
            unless defined $header[$i] && $header[$i] ne "";

        if (exists $idx{ $header[$i] }) {
            die "[ERROR] Duplicate column name found in input header: $header[$i]\n";
        }

        $idx{ $header[$i] } = $i;
    }

    return %idx;
}

sub row_hash {
    my ($row_ref, $header_ref, $idx_ref) = @_;

    my @row    = @$row_ref;
    my @header = @$header_ref;

    my %r;

    for my $h (@header) {
        my $i = $idx_ref->{$h};
        $r{$h} = defined $row[$i] && $row[$i] ne "" ? $row[$i] : "NA";
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

    for my $key (@keys) {
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

    my %seen_sample;
    my $first_sample = "";

    for my $row (@$rows_ref) {
        my %r = row_hash($row, $header_ref, $idx_ref);
        my $s = get_value(\%r, $sample_col);

        next if $s eq "NA" || $s eq "";

        $seen_sample{$s}++;

        if ($first_sample eq "") {
            $first_sample = $s;
        }
    }

    if ($first_sample eq "") {
        return "NA";
    }

    if (scalar(keys %seen_sample) > 1) {
        print STDERR "[WARN] Multiple sample names found in first column '$sample_col'.\n";
        print STDERR "[WARN] Output file prefix will use the first sample: $first_sample\n";
    }

    return $first_sample;
}

sub clean_tsv_value {
    my ($v) = @_;

    $v = "NA" unless defined $v;
    $v = "NA" if $v eq "";

    $v =~ s/\r/ /g;
    $v =~ s/\n/ /g;
    $v =~ s/\t/ /g;

    return $v;
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

            $key =~ s/^\s+|\s+$//g;
            $val =~ s/^\s+|\s+$//g;

            $val =~ s/^['"]//;
            $val =~ s/['"]$//;

            die "[ERROR] Empty config key found in $file\n"
                if $key eq "";

            $conf{$key} = $val;
        }
    }

    close $fh;

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
# Path helpers
# ============================================================

sub guess_project_root {
    my ($config_file) = @_;

    my $dir = dirname($config_file);

    if ($dir =~ m{/conf$}) {
        $dir =~ s{/conf$}{};
        return $dir;
    }

    return dirname($config_file);
}

sub resolve_path {
    my ($path, $project_root, $config_dir) = @_;

    die "[ERROR] Empty path\n"
        unless defined $path && $path ne "";

    if ($path =~ m{^/}) {
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

    die "[ERROR] Empty output directory\n"
        unless defined $dir && $dir ne "";

    if ($dir =~ m{^/}) {
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
           $hash{$b} <=> $hash{$a}
        || $a cmp $b
    } keys %hash;

    return $sorted[0];
}

sub unique {
    my @x = @_;

    my %seen;
    my @u;

    for my $v (@x) {
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

    return defined $x && $x =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
}

sub is_numeric_positive {
    my ($x) = @_;

    return 0 unless defined $x;
    return 0 if $x eq "" || $x eq "NA";
    return 0 unless $x =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;

    return $x > 0 ? 1 : 0;
}

sub usage {
    return <<'USAGE';
Usage:
  perl bin/annotate_candidates.pl \
    --config conf/hcm_exondel.example.conf \
    --input test_results/SAMPLE/04.candidates/SAMPLE.merged_candidates.tsv \
    --out test_results/SAMPLE/05.report/SAMPLE.annotated_candidates.tsv

Purpose:
  Annotate all candidates in merged_candidates.tsv using MANE RefSeq exon
  annotation.

Current merge logic:
  merge_evidence.pl uses split-read clusters as candidate events.
  Depth and discordant read-pair evidence are used as validation evidence.

This script:
  1. Does not filter by Evidence_Level.
  2. Does not split High/Other outputs.
  3. Annotates High, Moderate and Low candidates.
  4. Supports empty candidate input files with header only.
  5. Preserves key evidence fields from merge_evidence.pl.

Coordinate priority:
  1. Chrom + Best_Start + Best_End
  2. Chrom + Split_Start + Split_End
  3. Chrom + Merged_Start + Merged_End
  4. Chrom + Core_Start + Core_End

Config parameters used:
  REFSEQ_MANE_SELECT_EXON_TXT
  MIN_EXON_OVERLAP_FRACTION

Expected exon TXT format:
  Gene Transcript Exon Chrom Start End Strand
USAGE
}

