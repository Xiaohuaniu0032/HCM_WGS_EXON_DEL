#!/usr/bin/env perl
# -*- coding: utf-8 -*-

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path);

# ============================================================
# annotate_candidates.pl
#
# Annotate final candidate deletion regions with exon annotation.
#
# Input:
#   1. final_report.tsv from merge_evidence.pl
#   2. EXON_ANNOTATION_BED from config
#
# Output:
#   annotated_candidates.tsv
#
# Required command-line arguments:
#   --config
#   --input
#   --out
#
# Config parameters:
#   EXON_ANNOTATION_BED
#   MIN_EXON_OVERLAP_FRACTION
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

die "[ERROR] Config file not found: $config\n" unless -s $config;
die "[ERROR] Input final report not found: $input\n" unless -s $input;

my %CONF = read_config($config);

my $exon_bed = get_conf_required(\%CONF, "EXON_ANNOTATION_BED");
$exon_bed = abs_path($exon_bed);

die "[ERROR] EXON_ANNOTATION_BED file not found: $exon_bed\n" unless -s $exon_bed;

my $min_exon_overlap_fraction =
    get_conf_value(\%CONF, "MIN_EXON_OVERLAP_FRACTION", 0.20);

my $outdir = dirname($out);
make_path($outdir) unless -d $outdir;

print "[INFO] Candidate annotation started\n";
print "[INFO] Input final report          : $input\n";
print "[INFO] Exon annotation BED         : $exon_bed\n";
print "[INFO] Min exon overlap fraction   : $min_exon_overlap_fraction\n";
print "[INFO] Output                      : $out\n";

# ============================================================
# Step 1. Read exon annotation
# ============================================================

my @exons = read_exon_bed($exon_bed);

die "[ERROR] No valid exon records found in EXON_ANNOTATION_BED: $exon_bed\n"
    unless @exons;

# group by chr for speed
my %exons_by_chr;
foreach my $exon (@exons) {
    push @{ $exons_by_chr{ $exon->{chr} } }, $exon;
}

foreach my $chr (keys %exons_by_chr) {
    @{ $exons_by_chr{$chr} } = sort {
        $a->{start1} <=> $b->{start1}
        ||
        $a->{end1} <=> $b->{end1}
    } @{ $exons_by_chr{$chr} };
}

# ============================================================
# Step 2. Read final candidates and annotate
# ============================================================

my ($header_ref, $rows_ref) = read_tsv($input);
my @header = @$header_ref;
my @rows   = @$rows_ref;

my %idx = header_index(@header);

foreach my $required (qw(SampleID Chrom Start End Candidate_Length)) {
    die "[ERROR] Required column '$required' not found in input: $input\n"
        unless exists $idx{$required};
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

my $annotated_count = 0;
my $total_count = 0;

foreach my $row (@rows) {
    $total_count++;

    my %r = row_hash($row, \@header, \%idx);

    my $chr   = $r{Chrom};
    my $start = $r{Start};
    my $end   = $r{End};

    my $annotation = annotate_one_candidate(
        chr                       => $chr,
        start                     => $start,
        end                       => $end,
        exons_by_chr              => \%exons_by_chr,
        min_exon_overlap_fraction => $min_exon_overlap_fraction,
    );

    $annotated_count++ if $annotation->{annotation_status} eq "Annotated";

    my @original_values = map {
        defined $r{$_} ? $r{$_} : "NA"
    } @header;

    print $out_fh join("\t",
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

close $out_fh;

print "[INFO] Candidate annotation finished\n";
print "[INFO] Total candidates     : $total_count\n";
print "[INFO] Annotated candidates : $annotated_count\n";
print "[INFO] Output               : $out\n";

exit 0;


# ============================================================
# Annotation functions
# ============================================================

sub annotate_one_candidate {
    my %args = @_;

    my $chr   = $args{chr};
    my $start = $args{start};
    my $end   = $args{end};

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
        last if $exon->{start1} > $end;
        next if $exon->{end1} < $start;

        my $overlap_len = overlap_length($start, $end, $exon->{start1}, $exon->{end1});
        next if $overlap_len <= 0;

        my $exon_len = $exon->{end1} - $exon->{start1} + 1;
        my $frac = $overlap_len / $exon_len;

        next if $frac < $min_fraction;

        my $exon_label = join("|",
            $exon->{gene},
            $exon->{transcript},
            $exon->{exon}
        );

        push @overlapped_exons, $exon_label;

        $genes{ $exon->{gene} }++;
        $transcripts{ $exon->{transcript} }++;

        my $detail = join(":",
            $exon_label,
            $exon->{chr},
            $exon->{start1},
            $exon->{end1},
            sprintf("%.3f", $frac)
        );

        push @details, $detail;

        if ($start <= $exon->{start1} && $end >= $exon->{end1}) {
            push @fully_covered, $exon_label;
        }
        else {
            push @partial, $exon_label;
        }
    }

    if (!@overlapped_exons) {
        return empty_annotation("No exon overlapped");
    }

    my $gene = choose_most_frequent_key(%genes) || "NA";
    my $transcript = choose_most_frequent_key(%transcripts) || "NA";

    return {
        gene                         => $gene,
        transcript                   => $transcript,
        affected_exons               => join(",", unique(@overlapped_exons)),
        overlap_exon_count           => scalar(unique(@overlapped_exons)),
        fully_covered_exons          => @fully_covered ? join(",", unique(@fully_covered)) : "NA",
        partially_overlapped_exons   => @partial ? join(",", unique(@partial)) : "NA",
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
# Read exon BED
# ============================================================

sub read_exon_bed {
    my ($file) = @_;

    my @exons;

    open my $fh, "<", $file
        or die "[ERROR] Cannot open exon BED: $file\n";

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        my @f = split /\t/, $line;

        die "[ERROR] Invalid exon BED line: $line\n" unless @f >= 4;

        my ($chr, $start0, $end1, $name) = @f[0, 1, 2, 3];

        die "[ERROR] Invalid exon BED coordinate: $line\n"
            unless $start0 =~ /^\d+$/ && $end1 =~ /^\d+$/ && $end1 > $start0;

        my ($gene, $transcript, $exon) = parse_exon_name($name);

        push @exons, {
            chr        => $chr,
            start0     => $start0,
            start1     => $start0 + 1,
            end1       => $end1,
            name       => $name,
            gene       => $gene,
            transcript => $transcript,
            exon       => $exon,
        };
    }

    close $fh;

    return @exons;
}


sub parse_exon_name {
    my ($name) = @_;

    my ($gene, $transcript, $exon);

    if ($name =~ /\|/) {
        my @x = split /\|/, $name;
        $gene       = $x[0] || "NA";
        $transcript = $x[1] || "NA";
        $exon       = $x[2] || "NA";
    }
    elsif ($name =~ /:/) {
        my @x = split /:/, $name;
        $gene       = $x[0] || "NA";
        $transcript = $x[1] || "NA";
        $exon       = $x[2] || "NA";
    }
    else {
        $gene       = $name || "NA";
        $transcript = "NA";
        $exon       = "NA";
    }

    return ($gene, $transcript, $exon);
}


# ============================================================
# Generic TSV helpers
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


# ============================================================
# Misc helpers
# ============================================================

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


sub usage {
    return <<"USAGE";

Usage:

  perl bin/annotate_candidates.pl \\
    --config conf/hcm_exondel.conf \\
    --input results/SAMPLE001/05.report/SAMPLE001.final_report.tsv \\
    --out results/SAMPLE001/05.report/SAMPLE001.annotated_candidates.tsv

Required arguments:

  --config       Config file
  --input        Final report from merge_evidence.pl
  --out          Annotated candidate output

Config parameters used:

  EXON_ANNOTATION_BED
  MIN_EXON_OVERLAP_FRACTION

Expected EXON_ANNOTATION_BED format:

  chr18    34232239    34233000    FHOD3|NM_001281740.3|EX12

Output:

  annotated_candidates.tsv

USAGE
}

