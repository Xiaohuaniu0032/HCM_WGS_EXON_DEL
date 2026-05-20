#!/usr/bin/env perl
# -*- coding: utf-8 -*-

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path getcwd);

# ============================================================
# annotate_candidates.pl
#
# Annotate merged candidate deletion regions with MANE RefSeq
# exon annotation.
#
# Input:
#   merged_candidates.tsv or final_report.tsv from merge_evidence.pl
#
# Required input columns:
#   Chrom
#   Start
#   End
#
# Required when ANNOTATE_EVIDENCE_LEVEL is not All:
#   Evidence_Level
#
# Config parameters used:
#   REFSEQ_MANE_SELECT_EXON_TXT
#   MIN_EXON_OVERLAP_FRACTION
#   ANNOTATE_EVIDENCE_LEVEL
#
# Exon TXT format:
#   Gene    Transcript    Exon    Chrom    Start    End    Strand
#
# Coordinate system:
#   1-based closed interval
#
# Default behavior:
#   ANNOTATE_EVIDENCE_LEVEL=High
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

die "[ERROR] Config file not found: $config\n" unless defined $config && -s $config;
die "[ERROR] Input candidate file not found: $input\n" unless defined $input && -s $input;

my %CONF = read_config($config);

my $project_root = guess_project_root($config);
my $config_dir   = dirname($config);

my $exon_txt = get_conf_required(\%CONF, "REFSEQ_MANE_SELECT_EXON_TXT");
$exon_txt = resolve_path($exon_txt, $project_root, $config_dir);

die "[ERROR] REFSEQ_MANE_SELECT_EXON_TXT file not found: $exon_txt\n" unless -s $exon_txt;

my $min_exon_overlap_fraction = get_conf_value(
    \%CONF,
    "MIN_EXON_OVERLAP_FRACTION",
    0.20
);

die "[ERROR] MIN_EXON_OVERLAP_FRACTION must be numeric: $min_exon_overlap_fraction\n"
    unless is_number($min_exon_overlap_fraction);

die "[ERROR] MIN_EXON_OVERLAP_FRACTION must be > 0 and <= 1: $min_exon_overlap_fraction\n"
    unless $min_exon_overlap_fraction > 0 && $min_exon_overlap_fraction <= 1;

my $annotate_evidence_level = get_conf_value(
    \%CONF,
    "ANNOTATE_EVIDENCE_LEVEL",
    "High"
);

$annotate_evidence_level =~ s/^\s+//;
$annotate_evidence_level =~ s/\s+$//;

die "[ERROR] ANNOTATE_EVIDENCE_LEVEL is empty\n"
    if !defined $annotate_evidence_level || $annotate_evidence_level eq "";

my $annotate_all = 0;
if (lc($annotate_evidence_level) eq "all") {
    $annotate_all = 1;
}

my $outdir = dirname($out);
make_path($outdir) unless -d $outdir;

print "[INFO] Candidate annotation started\n";
print "[INFO] Config                    : $config\n";
print "[INFO] Project root              : $project_root\n";
print "[INFO] Input candidate file       : $input\n";
print "[INFO] Exon TXT                  : $exon_txt\n";
print "[INFO] Min exon overlap fraction : $min_exon_overlap_fraction\n";
print "[INFO] Annotate evidence level   : $annotate_evidence_level\n";
print "[INFO] Annotate all candidates   : $annotate_all\n";
print "[INFO] Output                    : $out\n";

# ============================================================
# Step 1. Read exon annotation
# ============================================================

my @exons = read_exon_txt($exon_txt);

die "[ERROR] No valid exon records found in REFSEQ_MANE_SELECT_EXON_TXT: $exon_txt\n"
    unless @exons;

my %exons_by_chr;

foreach my $exon (@exons) {
    push @{ $exons_by_chr{ $exon->{chrom} } }, $exon;
}

foreach my $chr (keys %exons_by_chr) {
    @{ $exons_by_chr{$chr} } = sort {
        $a->{start} <=> $b->{start} ||
        $a->{end}   <=> $b->{end}
    } @{ $exons_by_chr{$chr} };
}

# ============================================================
# Step 2. Read candidate file
# ============================================================

my ($header_ref, $rows_ref) = read_tsv($input);
my @header = @$header_ref;
my @rows   = @$rows_ref;
my %idx    = header_index(@header);

foreach my $required (qw(Chrom Start End)) {
    die "[ERROR] Required column '$required' not found in input: $input\n"
        unless exists $idx{$required};
}

if (!$annotate_all) {
    die "[ERROR] Required column 'Evidence_Level' not found in input: $input\n"
        unless exists $idx{Evidence_Level};
}

open my $out_fh, ">", $out
    or die "[ERROR] Cannot write output file: $out\n";

my @new_header = (
    @header,
    qw(
        Annotated_Gene
        Annotated_Transcript
        Affected_Exons
        Overlap_Exon_Count
        Fully_Covered_Exons
        Partially_Overlapped_Exons
        Exon_Overlap_Detail
        Annotation_Status
    )
);

print $out_fh join("\t", @new_header) . "\n";

my $total_count     = 0;
my $annotated_count = 0;
my $not_annotated   = 0;
my $skipped_count   = 0;

foreach my $row (@rows) {
    $total_count++;

    my %r = row_hash($row, \@header, \%idx);

    my $chr   = $r{Chrom};
    my $start = $r{Start};
    my $end   = $r{End};

    if (!defined $chr || $chr eq "" || $chr eq "NA") {
        die "[ERROR] Empty Chrom value at candidate row $total_count\n";
    }

    if (!is_integer($start) || !is_integer($end)) {
        die "[ERROR] Invalid Start/End at candidate row $total_count: Start=$start End=$end\n";
    }

    if ($start > $end) {
        ($start, $end) = ($end, $start);
    }

    my $annotation;

    if (!$annotate_all) {
        my $level = $r{Evidence_Level};

        if (!defined $level || $level ne $annotate_evidence_level) {
            $skipped_count++;

            $annotation = {
                gene                       => "NA",
                transcript                 => "NA",
                affected_exons             => "NA",
                overlap_exon_count         => 0,
                fully_covered_exons        => "NA",
                partially_overlapped_exons => "NA",
                exon_overlap_detail        => "Evidence_Level=$level; required=$annotate_evidence_level",
                annotation_status          => "Skipped_by_evidence_level_filter",
            };

            print_output_row($out_fh, \@header, \%r, $annotation);
            next;
        }
    }

    $annotation = annotate_one_candidate(
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

    print_output_row($out_fh, \@header, \%r, $annotation);
}

close $out_fh;

print "[INFO] Candidate annotation finished\n";
print "[INFO] Total candidates       : $total_count\n";
print "[INFO] Skipped candidates     : $skipped_count\n";
print "[INFO] Annotated candidates   : $annotated_count\n";
print "[INFO] Not annotated          : $not_annotated\n";
print "[INFO] Output                 : $out\n";

exit 0;

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

    my @header = split /\t/, $header;
    my %idx    = header_index(@header);

    foreach my $required (qw(Gene Transcript Exon Chrom Start End Strand)) {
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
# Output helper
# ============================================================

sub print_output_row {
    my ($fh, $header_ref, $row_hash_ref, $annotation) = @_;

    my @original_values = map {
        defined $row_hash_ref->{$_} ? $row_hash_ref->{$_} : "NA"
    } @$header_ref;

    print $fh join(
        "\t",
        @original_values,
        $annotation->{gene},
        $annotation->{transcript},
        $annotation->{affected_exons},
        $annotation->{overlap_exon_count},
        $annotation->{fully_covered_exons},
        $annotation->{partially_overlapped_exons},
        $annotation->{exon_overlap_detail},
        $annotation->{annotation_status},
    ) . "\n";
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

# ============================================================
# Misc helpers
# ============================================================

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
    --out    test_results/SAMPLE/04.candidates/SAMPLE.annotated_candidates.tsv

Required arguments:
  --config    Config file
  --input     Candidate TSV file from merge_evidence.pl
  --out       Annotated output TSV file

Config parameters used:
  REFSEQ_MANE_SELECT_EXON_TXT
  MIN_EXON_OVERLAP_FRACTION
  ANNOTATE_EVIDENCE_LEVEL

Recommended config:
  ANNOTATE_EVIDENCE_LEVEL=High

Supported ANNOTATE_EVIDENCE_LEVEL values:
  High
  Medium
  Low
  All

Default behavior:
  ANNOTATE_EVIDENCE_LEVEL=High

  This means only rows with Evidence_Level=High are annotated.
  Other rows are kept in output and marked as:
  Skipped_by_evidence_level_filter

If:
  ANNOTATE_EVIDENCE_LEVEL=All

  All candidates are annotated regardless of Evidence_Level.

Expected exon TXT format:
  Gene    Transcript    Exon    Chrom    Start    End    Strand

Example:
  perl bin/annotate_candidates.pl \\
    --config conf/hcm_exondel.example.conf \\
    --input test/test_results/25B09089386.final.merge/04.candidates/25B09089386.final.merge.merged_candidates.tsv \\
    --out test/test_results/25B09089386.final.merge/04.candidates/25B09089386.final.merge.annotated_candidates.tsv

USAGE
}

