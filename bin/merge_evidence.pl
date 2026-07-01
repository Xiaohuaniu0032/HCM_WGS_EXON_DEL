#!/usr/bin/env perl
# -*- coding: utf-8 -*-
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename qw(dirname);
use File::Path qw(make_path);

# ============================================================
# merge_evidence.pl
#
# Function:
#   Split-read-centered evidence merging for exon-level deletion
#   candidates.
#
# Core logic:
#   For each row in split_reads.clusters.tsv:
#
#     1. Use Cluster_Start / Cluster_End as the split DEL interval.
#
#     2. Check depth_candidates.cusum.all.tsv.
#        If the union of overlapping depth intervals covers at least
#        --min-overlap fraction of the split DEL bases, this split DEL
#        is considered supported by depth evidence.
#
#     3. Check discordant_reads.tsv.
#        If any single discordant interval covers at least --min-overlap
#        fraction of the split DEL bases, this split DEL is considered
#        supported by discordant evidence.
#
# Evidence level:
#   High     = Split + Depth + Discordant
#   Moderate = Split + Depth
#   Low      = Split without Depth
#
# Notes:
#   1. Coordinates are treated as 1-based closed intervals.
#   2. Depth evidence should be CUSUM interval-level evidence:
#        SAMPLE.depth_candidates.cusum.all.tsv
#   3. This script does not read config directly.
#      HCMExonDel.pl reads EVIDENCE_OVERLAP_FRACTION from config
#      and passes it to this script as --min-overlap.
# ============================================================

my %opt;

GetOptions(
    "sample=s"      => \$opt{sample},
    "depth=s"       => \$opt{depth},
    "split=s"       => \$opt{split},
    "discordant=s"  => \$opt{discordant},
    "min-overlap=f" => \$opt{min_overlap},
    "out=s"         => \$opt{out},
    "help|h"        => \$opt{help},
) or die usage();

if ($opt{help}) {
    print usage();
    exit 0;
}

for my $k (qw/sample depth split discordant min_overlap out/) {
    die "[ERROR] Missing required option --" . option_name($k) . "\n" . usage()
        unless defined $opt{$k} && $opt{$k} ne "";
}

validate_fraction("min-overlap", $opt{min_overlap});

die "[ERROR] Split cluster file not found: $opt{split}\n"
    unless -e $opt{split};

die "[ERROR] Depth evidence file not found: $opt{depth}\n"
    unless -e $opt{depth};

die "[ERROR] Discordant evidence file not found: $opt{discordant}\n"
    unless -e $opt{discordant};

log_msg("Merge evidence started");
log_msg("Sample               : $opt{sample}");
log_msg("Split evidence file  : $opt{split}");
log_msg("Depth evidence file  : $opt{depth}");
log_msg("Discordant file      : $opt{discordant}");
log_msg("Min overlap fraction : $opt{min_overlap}");
log_msg("Output               : $opt{out}");

my @split_records      = read_split_records($opt{split}, $opt{sample});
my @depth_records      = read_interval_records($opt{depth}, $opt{sample}, "Depth");
my @discordant_records = read_interval_records($opt{discordant}, $opt{sample}, "Discordant");

log_msg("Split records loaded      : " . scalar(@split_records));
log_msg("Depth records loaded      : " . scalar(@depth_records));
log_msg("Discordant records loaded : " . scalar(@discordant_records));

my %depth_by_chrom      = index_by_chrom(@depth_records);
my %discordant_by_chrom = index_by_chrom(@discordant_records);

# If no split records are present, write header-only output.
if (!@split_records) {
    write_output($opt{out}, []);
    log_msg("No split records found. Header-only output written.");
    log_msg("Merge evidence finished");
    exit 0;
}

my @out_records;

for my $sp (@split_records) {
    my $chrom = $sp->{Chrom};
    my $s     = $sp->{Start};
    my $e     = $sp->{End};
    my $len   = interval_len($s, $e);

    # ========================================================
    # Depth support
    # ========================================================
    # Inclusion condition:
    #   same chromosome and overlap_bases > 0
    #
    # Support condition:
    #   union of all clipped overlapping depth intervals covers
    #   >= min-overlap fraction of split DEL bases.
    # ========================================================

    my @depth_overlaps;

    for my $dp (@{ $depth_by_chrom{$chrom} || [] }) {
        my $ov = interval_overlap_len($s, $e, $dp->{Start}, $dp->{End});
        next unless $ov > 0;

        push @depth_overlaps, $dp;
    }

    my ($depth_cov_bases, $depth_cov_frac, $depth_union_ranges) =
        covered_bases_by_union(\@depth_overlaps, $chrom, $s, $e);

    my $depth_support = ($depth_cov_frac >= $opt{min_overlap}) ? 1 : 0;

    my ($depth_overlap_ranges, $depth_overlap_bases_by_record) =
        clipped_overlap_detail(\@depth_overlaps, $s, $e);

    # ========================================================
    # Discordant support
    # ========================================================
    # Inclusion condition:
    #   same chromosome and overlap_bases > 0
    #
    # Support condition:
    #   at least one single discordant interval covers
    #   >= min-overlap fraction of split DEL bases.
    # ========================================================

    my @discordant_hits;

    my $best_discordant_bases = 0;
    my $best_discordant_frac  = 0;
    my $best_discordant_range = "NA";
    my $best_discordant_id    = "NA";

    for my $dr (@{ $discordant_by_chrom{$chrom} || [] }) {
        my $ov = interval_overlap_len($s, $e, $dr->{Start}, $dr->{End});
        next unless $ov > 0;

        my $fr = $len > 0 ? $ov / $len : 0;

        if ($fr > $best_discordant_frac) {
            $best_discordant_frac  = $fr;
            $best_discordant_bases = $ov;
            $best_discordant_range = $dr->{Chrom} . ":" . $dr->{Start} . "-" . $dr->{End};
            $best_discordant_id    = $dr->{Source_ID};
        }

        if ($fr >= $opt{min_overlap}) {
            push @discordant_hits, $dr;
        }
    }

    my $discordant_support = @discordant_hits ? 1 : 0;

    my ($discordant_overlap_ranges, $discordant_overlap_bases_by_record) =
        clipped_overlap_detail(\@discordant_hits, $s, $e);

    # ========================================================
    # Evidence level
    # ========================================================

    my $level;

    if ($depth_support && $discordant_support) {
        $level = "High";
    }
    elsif ($depth_support) {
        $level = "Moderate";
    }
    else {
        $level = "Low";
    }

    my @evidence_types = ("Split");
    push @evidence_types, "Depth"      if $depth_support;
    push @evidence_types, "Discordant" if $discordant_support;

    my $comment = make_comment(
        level                => $level,
        min_overlap          => $opt{min_overlap},
        depth_cov_frac       => $depth_cov_frac,
        best_discordant_frac => $best_discordant_frac,
    );

    push @out_records, {
        SampleID                              => $sp->{SampleID},
        Gene                                  => $sp->{Gene},
        Chrom                                 => $chrom,
        Cluster_ID                            => $sp->{Source_ID},

        Split_Start                           => $s,
        Split_End                             => $e,
        Split_Size                            => $len,
        Split_Reads                           => meta_value_any($sp, qw/Support_Reads Split_Reads Read_Count N_Reads/),
        Split_Records                         => meta_value_any($sp, qw/Support_Records Split_Records Record_Count N_Records/),
        Split_Status                          => meta_value_any($sp, qw/Cluster_Status Status Candidate_Status/),

        Depth_Support                         => $depth_support ? "Yes" : "No",
        Depth_Covered_Bases                   => $depth_cov_bases,
        Depth_Coverage_Fraction               => sprintf("%.4f", $depth_cov_frac),
        Depth_Record_Count                    => scalar(@depth_overlaps),
        Depth_Range                           => ranges_for_records(@depth_overlaps),
        Depth_Overlap_Range                   => $depth_overlap_ranges,
        Depth_Overlap_Bases_By_Record         => $depth_overlap_bases_by_record,
        Depth_Union_Range                     => $depth_union_ranges,

        Discordant_Read_Support               => $discordant_support ? "Yes" : "No",
        Discordant_Overlap_Bases              => $best_discordant_bases,
        Discordant_Overlap_Fraction           => sprintf("%.4f", $best_discordant_frac),
        Discordant_Record_Count               => scalar(@discordant_hits),
        Discordant_Range                      => ranges_for_records(@discordant_hits),
        Discordant_Overlap_Range              => $discordant_overlap_ranges,
        Discordant_Overlap_Bases_By_Record    => $discordant_overlap_bases_by_record,
        Discordant_Cluster_IDs                => uniq_values(map { $_->{Source_ID} } @discordant_hits),
        Discordant_Reads                      => sum_meta_numeric_any(
                                                     \@discordant_hits,
                                                     qw/Discordant_Reads Support_Reads Read_Count N_Reads/
                                                 ),
        Median_Insert_Size                    => uniq_values(
                                                     map {
                                                         meta_value_any(
                                                             $_,
                                                             qw/Median_Insert_Size Median_TLEN Median_InsertSize Insert_Size/
                                                         )
                                                     } @discordant_hits
                                                 ),
        Best_Discordant_Cluster_ID            => $best_discordant_id,
        Best_Discordant_Range                 => $best_discordant_range,

        Min_Overlap_Fraction                  => $opt{min_overlap},
        Evidence_Count                        => scalar(@evidence_types),
        Evidence_Level                        => $level,
        Evidence_Types                        => join(",", @evidence_types),

        Best_Start                            => $s,
        Best_End                              => $e,
        Best_Size                             => $len,
        Best_Evidence                         => "Split",

        Candidate_Status                      => candidate_status($level),
        Source_Records                        => join_source_records($sp, \@depth_overlaps, \@discordant_hits),
        Comment                               => $comment,
    };
}

@out_records = sort {
       chrom_cmp($a->{Chrom}, $b->{Chrom})
    || $a->{Split_Start} <=> $b->{Split_Start}
    || $a->{Split_End}   <=> $b->{Split_End}
    || $a->{Cluster_ID}  cmp $b->{Cluster_ID}
} @out_records;

write_output($opt{out}, \@out_records);

my $n_high = scalar(grep { $_->{Evidence_Level} eq "High" } @out_records);
my $n_mod  = scalar(grep { $_->{Evidence_Level} eq "Moderate" } @out_records);
my $n_low  = scalar(grep { $_->{Evidence_Level} eq "Low" } @out_records);

log_msg("Output candidates written : " . scalar(@out_records));
log_msg("High=$n_high Moderate=$n_mod Low=$n_low");
log_msg("Merge evidence finished");

exit 0;

# ============================================================
# Readers
# ============================================================

sub read_split_records {
    my ($file, $sample_filter) = @_;

    return () unless -s $file;

    my ($cols, $rows) = read_table($file);
    my %idx = index_cols($cols);

    my $sample_col = find_col(\%idx, qw/Sample SampleID/);
    my $id_col     = require_any_col($file, \%idx, qw/Cluster_ID ID/);
    my $chrom_col  = require_any_col($file, \%idx, qw/SV_Chrom Chrom CHROM/);
    my $start_col  = require_any_col($file, \%idx, qw/Cluster_Start Start START/);
    my $end_col    = require_any_col($file, \%idx, qw/Cluster_End End END/);
    my $gene_col   = find_col(\%idx, qw/Scan_Genes Gene Genes Target_Gene Target_Name/);

    my @out;

    for my $row (@$rows) {
        my $sample = defined $sample_col
            ? normalize_value($row->{$sample_col})
            : $sample_filter;

        next if defined $sample_filter
             && $sample_filter ne ""
             && $sample ne $sample_filter;

        my ($s, $e) = normalize_interval($row->{$start_col}, $row->{$end_col}, $file);

        my $gene = defined $gene_col
            ? normalize_value($row->{$gene_col})
            : "NA";

        push @out, {
            Type      => "Split",
            SampleID  => $sample,
            Gene      => $gene,
            Chrom     => normalize_value($row->{$chrom_col}),
            Start     => $s,
            End       => $e,
            Source_ID => normalize_value($row->{$id_col}),
            Meta      => { %$row },
        };
    }

    return @out;
}

sub read_interval_records {
    my ($file, $sample_filter, $type) = @_;

    return () unless -s $file;

    my ($cols, $rows) = read_table($file);
    my %idx = index_cols($cols);

    my $sample_col     = find_col(\%idx, qw/SampleID Sample/);
    my $id_col         = find_col(\%idx, qw/Cluster_ID Del_Window_IDs Depth_File ID/);
    my $chrom_col      = require_any_col($file, \%idx, qw/Chrom CHROM SV_Chrom/);
    my $start_col      = require_any_col($file, \%idx, qw/Start START Cluster_Start/);
    my $end_col        = require_any_col($file, \%idx, qw/End END Cluster_End/);
    my $gene_col       = find_col(\%idx, qw/Gene Scan_Genes Genes Target_Gene Target_Name/);
    my $type_col       = find_col(\%idx, qw/Type TYPE SV_Type SVTYPE/);
    my $depth_file_col = find_col(\%idx, qw/Depth_File DepthFile File/);

    my @out;
    my $n = 0;

    for my $row (@$rows) {
        $n++;

        if (defined $sample_col) {
            my $sample = normalize_value($row->{$sample_col});

            next if defined $sample_filter
                 && $sample_filter ne ""
                 && $sample ne $sample_filter;
        }

        if ($type eq "Depth" && defined $type_col) {
            my $svtype = normalize_value($row->{$type_col});
            next if $svtype ne "NA" && $svtype !~ /DEL/i;
        }

        my ($s, $e) = normalize_interval($row->{$start_col}, $row->{$end_col}, $file);

        my $gene = defined $gene_col
            ? normalize_value($row->{$gene_col})
            : "NA";

        if ($gene eq "NA" && defined $depth_file_col) {
            $gene = guess_gene_from_depth_file($row->{$depth_file_col});
        }

        my $id = defined $id_col ? normalize_value($row->{$id_col}) : "NA";
        $id = $type . "_" . $n if $id eq "NA";

        push @out, {
            Type      => $type,
            SampleID  => $sample_filter,
            Gene      => $gene,
            Chrom     => normalize_value($row->{$chrom_col}),
            Start     => $s,
            End       => $e,
            Source_ID => $id,
            Meta      => { %$row },
        };
    }

    return @out;
}

sub read_table {
    my ($file) = @_;

    open my $fh, "<", $file
        or die "[ERROR] Cannot open $file: $!\n";

    my $header = <$fh>;

    die "[ERROR] Empty file: $file\n"
        unless defined $header;

    chomp $header;
    $header =~ s/\r$//;

    my @cols = split_line($header);

    die "[ERROR] No header columns found in $file\n"
        unless @cols;

    my @rows;

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;

        my @v = split_line($line, scalar(@cols));

        my %row;

        for my $i (0 .. $#cols) {
            $row{$cols[$i]} = defined $v[$i] && $v[$i] ne "" ? $v[$i] : "NA";
        }

        push @rows, \%row;
    }

    close $fh;

    return (\@cols, \@rows);
}

sub split_line {
    my ($line, $limit) = @_;

    if ($line =~ /\t/) {
        return split /\t/, $line, defined $limit ? $limit : 0;
    }

    return split /\s+/, $line, defined $limit ? $limit : 0;
}

# ============================================================
# Evidence logic
# ============================================================

sub covered_bases_by_union {
    my ($records_ref, $chrom, $query_start, $query_end) = @_;

    my @iv;

    for my $r (@$records_ref) {
        my $s = max_num($query_start, $r->{Start});
        my $e = min_num($query_end,   $r->{End});

        next if $s > $e;

        push @iv, [$s, $e];
    }

    return (0, 0, "NA") unless @iv;

    @iv = sort {
           $a->[0] <=> $b->[0]
        || $a->[1] <=> $b->[1]
    } @iv;

    my @merged;
    my ($cs, $ce) = @{ shift @iv };

    for my $x (@iv) {
        my ($s, $e) = @$x;

        if ($s <= $ce + 1) {
            $ce = $e if $e > $ce;
        }
        else {
            push @merged, [$cs, $ce];
            ($cs, $ce) = ($s, $e);
        }
    }

    push @merged, [$cs, $ce];

    my $covered = 0;

    for my $m (@merged) {
        $covered += $m->[1] - $m->[0] + 1;
    }

    my $query_len = interval_len($query_start, $query_end);
    my $frac = $query_len > 0 ? $covered / $query_len : 0;

    my $ranges = join(
        ",",
        map { $chrom . ":" . $_->[0] . "-" . $_->[1] } @merged
    );

    return ($covered, $frac, $ranges);
}

sub clipped_overlap_detail {
    my ($records_ref, $query_start, $query_end) = @_;

    my @ranges;
    my @bases;

    for my $r (@$records_ref) {
        my $s = max_num($query_start, $r->{Start});
        my $e = min_num($query_end,   $r->{End});

        next if $s > $e;

        my $b = $e - $s + 1;

        push @ranges, $r->{Chrom} . ":" . $s . "-" . $e;
        push @bases,  $b;
    }

    my $range_text = @ranges ? join(",", @ranges) : "NA";
    my $base_text  = @bases  ? join(",", @bases)  : "NA";

    return ($range_text, $base_text);
}

sub interval_overlap_len {
    my ($s1, $e1, $s2, $e2) = @_;

    my $s = max_num($s1, $s2);
    my $e = min_num($e1, $e2);

    return 0 if $s > $e;

    return $e - $s + 1;
}

sub interval_len {
    my ($s, $e) = @_;

    return $e >= $s ? $e - $s + 1 : 0;
}

# ============================================================
# Output helpers
# ============================================================

sub write_output {
    my ($file, $records_ref) = @_;

    my $dir = dirname($file);
    make_path($dir) if defined $dir && $dir ne "." && !-d $dir;

    my @header = qw/
        SampleID
        Gene
        Chrom
        Cluster_ID

        Split_Start
        Split_End
        Split_Size
        Split_Reads
        Split_Records
        Split_Status

        Depth_Support
        Depth_Covered_Bases
        Depth_Coverage_Fraction
        Depth_Record_Count
        Depth_Range
        Depth_Overlap_Range
        Depth_Overlap_Bases_By_Record
        Depth_Union_Range

        Discordant_Read_Support
        Discordant_Overlap_Bases
        Discordant_Overlap_Fraction
        Discordant_Record_Count
        Discordant_Range
        Discordant_Overlap_Range
        Discordant_Overlap_Bases_By_Record
        Discordant_Cluster_IDs
        Discordant_Reads
        Median_Insert_Size
        Best_Discordant_Cluster_ID
        Best_Discordant_Range

        Min_Overlap_Fraction
        Evidence_Count
        Evidence_Level
        Evidence_Types

        Best_Start
        Best_End
        Best_Size
        Best_Evidence

        Candidate_Status
        Source_Records
        Comment
    /;

    open my $out, ">", $file
        or die "[ERROR] Cannot write $file: $!\n";

    print $out join("\t", @header), "\n";

    for my $r (@$records_ref) {
        my @v = map {
            defined $r->{$_} && $r->{$_} ne "" ? $r->{$_} : "NA"
        } @header;

        print $out join("\t", @v), "\n";
    }

    close $out;
}

sub make_comment {
    my %arg = @_;

    my $level                = $arg{level};
    my $min_overlap          = $arg{min_overlap};
    my $depth_cov_frac       = sprintf("%.4f", $arg{depth_cov_frac});
    my $best_discordant_frac = sprintf("%.4f", $arg{best_discordant_frac});

    if ($level eq "High") {
        return "Split DEL is supported by both depth and discordant evidence; depth_fraction=$depth_cov_frac; best_discordant_fraction=$best_discordant_frac";
    }

    if ($level eq "Moderate") {
        return "Split DEL is supported by depth evidence only; depth_fraction=$depth_cov_frac; best_discordant_fraction=$best_discordant_frac";
    }

    return "Split DEL lacks sufficient depth support; depth_fraction=$depth_cov_frac; required_fraction=$min_overlap";
}

sub candidate_status {
    my ($level) = @_;

    return "High_confidence_candidate"     if $level eq "High";
    return "Moderate_confidence_candidate" if $level eq "Moderate";
    return "Low_confidence_candidate";
}

sub ranges_for_records {
    my @records = @_;

    return "NA" unless @records;

    return uniq_values(
        map {
            $_->{Chrom} . ":" . $_->{Start} . "-" . $_->{End}
        } @records
    );
}

sub join_source_records {
    my ($split, $depth_ref, $discordant_ref) = @_;

    my @x;

    push @x, source_label($split);
    push @x, map { source_label($_) } @$depth_ref;
    push @x, map { source_label($_) } @$discordant_ref;

    return uniq_values(@x);
}

sub source_label {
    my ($r) = @_;

    return join(
        ":",
        $r->{Type},
        $r->{Source_ID},
        $r->{Chrom} . ":" . $r->{Start} . "-" . $r->{End}
    );
}

sub sum_meta_numeric_any {
    my ($records_ref, @fields) = @_;

    my $sum = 0;
    my $has = 0;
    my %seen;

    for my $r (@$records_ref) {
        next if $seen{$r->{Source_ID}}++;

        my $v = meta_value_any($r, @fields);

        next unless defined $v && $v =~ /^-?\d+(?:\.\d+)?$/;

        $sum += $v;
        $has = 1;
    }

    return $has ? $sum : "NA";
}

sub meta_value_any {
    my ($r, @fields) = @_;

    return "NA" unless defined $r && exists $r->{Meta};

    for my $field (@fields) {
        return normalize_value($r->{Meta}{$field})
            if exists $r->{Meta}{$field};

        my $lc = lc $field;

        for my $k (keys %{ $r->{Meta} }) {
            return normalize_value($r->{Meta}{$k})
                if lc($k) eq $lc;
        }
    }

    return "NA";
}

# ============================================================
# General helpers
# ============================================================

sub index_by_chrom {
    my @records = @_;

    my %h;

    for my $r (@records) {
        push @{ $h{ $r->{Chrom} } }, $r;
    }

    return %h;
}

sub index_cols {
    my ($cols_ref) = @_;

    my %idx;

    for my $c (@$cols_ref) {
        $idx{$c}    = $c;
        $idx{lc $c} = $c;
    }

    return %idx;
}

sub find_col {
    my ($idx_ref, @candidates) = @_;

    for my $c (@candidates) {
        return $idx_ref->{$c}    if exists $idx_ref->{$c};
        return $idx_ref->{lc $c} if exists $idx_ref->{lc $c};
    }

    return undef;
}

sub require_any_col {
    my ($file, $idx_ref, @candidates) = @_;

    my $c = find_col($idx_ref, @candidates);

    die "[ERROR] Required column not found in $file. Accepted names: "
      . join(",", @candidates) . "\n"
        unless defined $c;

    return $c;
}

sub normalize_interval {
    my ($s, $e, $file) = @_;

    $s = normalize_value($s);
    $e = normalize_value($e);

    die "[ERROR] Invalid interval in $file: Start=$s End=$e\n"
        unless $s =~ /^\d+$/ && $e =~ /^\d+$/;

    ($s, $e) = ($e, $s) if $s > $e;

    return ($s + 0, $e + 0);
}

sub normalize_value {
    my ($v) = @_;

    return "NA" unless defined $v;

    $v =~ s/^\s+//;
    $v =~ s/\s+$//;

    return $v eq "" ? "NA" : $v;
}

sub guess_gene_from_depth_file {
    my ($v) = @_;

    $v = normalize_value($v);

    return "NA" if $v eq "NA";

    my $base = $v;
    $base =~ s{.*/}{};
    $base =~ s/\.depth\.tsv$//;

    # Example:
    # HCM_FHOD3_DEL_1000_2000_rep2_FHOD3_NM_001281740_chr18_...
    if ($base =~ /_([A-Za-z0-9.-]+)_N[MR]_\d+/) {
        return $1;
    }

    if ($base =~ /^HCM_([A-Za-z0-9.-]+)_DEL/) {
        return $1;
    }

    return "NA";
}

sub uniq_values {
    my @v = @_;

    my %seen;
    my @out;

    for my $x (@v) {
        $x = normalize_value($x);

        next if $x eq "NA";
        next if $seen{$x}++;

        push @out, $x;
    }

    return @out ? join(",", @out) : "NA";
}

sub chrom_cmp {
    my ($a, $b) = @_;

    my $ra = chrom_rank($a);
    my $rb = chrom_rank($b);

    return $ra <=> $rb || $a cmp $b;
}

sub chrom_rank {
    my ($c) = @_;

    $c = normalize_value($c);
    $c =~ s/^chr//i;

    return $c if $c =~ /^\d+$/;
    return 23 if uc($c) eq "X";
    return 24 if uc($c) eq "Y";
    return 25 if uc($c) =~ /^(M|MT)$/;

    return 1000;
}

sub min_num {
    return $_[0] < $_[1] ? $_[0] : $_[1];
}

sub max_num {
    return $_[0] > $_[1] ? $_[0] : $_[1];
}

sub validate_fraction {
    my ($name, $v) = @_;

    die "[ERROR] $name must be numeric: $v\n"
        unless defined $v
        && $v =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;

    die "[ERROR] $name must be > 0 and <= 1: $v\n"
        unless $v > 0 && $v <= 1;

    return 1;
}

sub option_name {
    my ($k) = @_;

    $k =~ s/_/-/g;

    return $k;
}

sub log_msg {
    my ($msg) = @_;

    print STDERR "[INFO] $msg\n";
}

sub usage {
    return <<'USAGE';
Usage:
  perl merge_evidence.pl \
    --sample SAMPLE_ID \
    --depth SAMPLE.depth_candidates.cusum.all.tsv \
    --split SAMPLE.split_reads.clusters.tsv \
    --discordant SAMPLE.discordant_reads.tsv \
    --min-overlap 0.90 \
    --out SAMPLE.merged_candidates.tsv

Required:
  --sample
      Sample ID.

  --depth
      CUSUM depth evidence file:
        SAMPLE.depth_candidates.cusum.all.tsv

  --split
      Split-read cluster file:
        SAMPLE.split_reads.clusters.tsv

  --discordant
      Discordant-read evidence file:
        SAMPLE.discordant_reads.tsv

  --min-overlap
      Minimum fraction of split DEL bases covered by depth or discordant
      evidence.

      In the HCMExonDel pipeline, this value is read from config:
        EVIDENCE_OVERLAP_FRACTION

      and then passed to this script by HCMExonDel.pl.

  --out
      Output merged candidate TSV.

Core logic:
  For each split DEL:
    1. Check whether depth evidence covers enough split DEL bases.
       Multiple overlapping depth intervals are merged by union first.

    2. Check whether any single discordant interval covers enough
       split DEL bases.

    3. Evidence level:
         High     = Split + Depth + Discordant
         Moderate = Split + Depth
         Low      = Split without Depth

Coordinates:
  1-based closed intervals.

Example:
  perl bin/merge_evidence.pl \
    --sample 25B09089386 \
    --depth test/test_results/25B09089386/01.depth/25B09089386.depth_candidates.cusum.all.tsv \
    --split test/test_results/25B09089386/02.split_reads/25B09089386.split_reads.clusters.tsv \
    --discordant test/test_results/25B09089386/03.discordant_reads/25B09089386.discordant_reads.tsv \
    --min-overlap 0.90 \
    --out test/test_results/25B09089386/04.candidates/25B09089386.merged_candidates.tsv
USAGE
}

