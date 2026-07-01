#!/usr/bin/env perl
# -*- coding: utf-8 -*-

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename qw(dirname);
use File::Path qw(make_path);

# ------------------------------------------------------------
# merge_evidence.pl
#
# Split-read-centered evidence merging logic.
#
# Core rule:
#   For each passed split-read cluster, evaluate whether the split interval
#   is supported by:
#     1. CUSUM depth evidence
#     2. Discordant read-pair evidence
#
# Confidence level:
#   High     = Split + Depth + Discordant
#   Moderate = Split + Depth only
#   Low      = Split without Depth, regardless of discordant support
#
# Important:
#   Split-read cluster is the candidate backbone.
#   Depth and discordant signals are validation evidence.
#
# Empty-result behavior:
#   If no passed split-read cluster exists, this script writes a header-only
#   merged_candidates.tsv and exits normally.
#
# Coordinates:
#   All coordinates are treated as 1-based closed intervals.
# ------------------------------------------------------------

my %opt = (
    min_overlap => 0.90,
);

GetOptions(
    'sample=s'      => \$opt{sample},
    'depth=s'       => \$opt{depth},
    'split=s'       => \$opt{split},
    'discordant=s'  => \$opt{discordant},
    'out=s'         => \$opt{out},
    'min-overlap=f' => \$opt{min_overlap},
    'help'          => \$opt{help},
) or die usage();

if ($opt{help}) {
    print usage();
    exit 0;
}

for my $k (qw/sample depth split discordant out/) {
    die "[ERROR] Missing required option --$k\n" . usage()
        unless defined $opt{$k} && $opt{$k} ne '';
}

check_existing_file($opt{depth},      'Depth file');
check_existing_file($opt{split},      'Split cluster file');
check_existing_file($opt{discordant}, 'Discordant file');

die "[ERROR] --min-overlap must be >0 and <=1\n"
    unless defined $opt{min_overlap}
        && $opt{min_overlap} > 0
        && $opt{min_overlap} <= 1;

log_msg("Merge evidence started");
log_msg("Sample     : $opt{sample}");
log_msg("Depth file : $opt{depth}");
log_msg("Split file : $opt{split}");
log_msg("Discordant : $opt{discordant}");
log_msg("Min overlap: $opt{min_overlap}");
log_msg("Output     : $opt{out}");

my @split_records      = read_split_records($opt{split}, $opt{sample});
my @depth_records      = read_interval_records($opt{depth}, $opt{sample}, 'Depth');
my @discordant_records = read_interval_records($opt{discordant}, $opt{sample}, 'Discordant');

log_msg("Split clusters loaded : " . scalar(@split_records));
log_msg("Depth records loaded   : " . scalar(@depth_records));
log_msg("Discordant loaded      : " . scalar(@discordant_records));

if (!@split_records) {
    log_msg("No passed split-read cluster found. Write header-only output and exit normally.");
    write_output($opt{out}, []);
    log_msg("Output candidates written: 0");
    log_msg("High=0 Moderate=0 Low=0");
    log_msg("Merge evidence finished");
    exit 0;
}

my %depth_by_chrom      = index_by_chrom(@depth_records);
my %discordant_by_chrom = index_by_chrom(@discordant_records);

my @out_records;

for my $sp (@split_records) {
    my $chrom = $sp->{Chrom};
    my $s     = $sp->{Start};
    my $e     = $sp->{End};
    my $len   = interval_len($s, $e);

    die "[ERROR] Invalid split interval length for "
      . "$chrom:$s-$e in split record $sp->{Source_ID}\n"
        unless $len > 0;

    # ------------------------------
    # Depth support
    # ------------------------------
    my @depth_candidates = records_for_chrom(\%depth_by_chrom, $chrom);

    my @depth_overlaps = grep {
        interval_overlap_len($s, $e, $_->{Start}, $_->{End}) > 0
    } @depth_candidates;

    my ($depth_cov_bases, $depth_cov_frac) = covered_bases_by_union(
        \@depth_overlaps,
        $s,
        $e
    );

    my $depth_support = ($depth_cov_frac >= $opt{min_overlap}) ? 1 : 0;

    # ------------------------------
    # Discordant support
    # ------------------------------
    # Conservative rule:
    #   Require at least one single discordant-read cluster to cover
    #   >= min_overlap of the split interval.
    my @discordant_candidates = records_for_chrom(\%discordant_by_chrom, $chrom);

    my @discordant_overlaps = grep {
        interval_overlap_len($s, $e, $_->{Start}, $_->{End}) > 0
    } @discordant_candidates;

    my @discordant_hits;
    my $best_discordant_bases = 0;
    my $best_discordant_frac  = 0;

    for my $dr (@discordant_overlaps) {
        my $ov = interval_overlap_len($s, $e, $dr->{Start}, $dr->{End});
        my $fr = $len > 0 ? $ov / $len : 0;

        if ($fr > $best_discordant_frac) {
            $best_discordant_frac  = $fr;
            $best_discordant_bases = $ov;
        }

        push @discordant_hits, $dr if $fr >= $opt{min_overlap};
    }

    my $discordant_support = @discordant_hits ? 1 : 0;

    # ------------------------------
    # Evidence level
    # ------------------------------
    my $level;

    if ($depth_support && $discordant_support) {
        $level = 'High';
    }
    elsif ($depth_support) {
        $level = 'Moderate';
    }
    else {
        $level = 'Low';
    }

    my @evidence_types = ('Split');

    push @evidence_types, 'Depth'      if $depth_support;
    push @evidence_types, 'Discordant' if $discordant_support;

    my $comment = make_comment(
        $level,
        $depth_cov_frac,
        $best_discordant_frac,
        $opt{min_overlap}
    );

    push @out_records, {
        SampleID                     => $sp->{SampleID},
        Gene                         => $sp->{Gene},
        Chrom                        => $chrom,
        Cluster_ID                   => $sp->{Source_ID},

        Split_Start                  => $s,
        Split_End                    => $e,
        Split_Size                   => $len,
        Split_Reads                  => meta_value($sp, 'Support_Reads'),
        Split_Records                => meta_value($sp, 'Support_Records'),

        Depth_Support                => $depth_support ? 'Yes' : 'No',
        Depth_Covered_Bases          => $depth_cov_bases,
        Depth_Coverage_Fraction      => sprintf('%.4f', $depth_cov_frac),
        Depth_Record_Count           => scalar(@depth_overlaps),
        Depth_Range                  => ranges_for_records(@depth_overlaps),

        Discordant_Read_Support      => $discordant_support ? 'Yes' : 'No',
        Discordant_Overlap_Bases     => $best_discordant_bases,
        Discordant_Overlap_Fraction  => sprintf('%.4f', $best_discordant_frac),
        Discordant_Record_Count      => scalar(@discordant_hits),
        Discordant_Range             => ranges_for_records(@discordant_hits),
        Discordant_Cluster_IDs       => uniq_values(map { $_->{Source_ID} } @discordant_hits),
        Discordant_Reads             => sum_meta_numeric(\@discordant_hits, 'Discordant_Reads'),
        Median_Insert_Size           => uniq_values(map { meta_value($_, 'Median_Insert_Size') } @discordant_hits),

        Evidence_Count               => scalar(@evidence_types),
        Evidence_Level               => $level,
        Evidence_Types               => join(',', @evidence_types),

        Best_Start                   => $s,
        Best_End                     => $e,
        Best_Size                    => $len,
        Best_Evidence                => 'Split',

        Candidate_Status             => candidate_status($level),
        Source_Records               => join_source_records($sp, \@depth_overlaps, \@discordant_hits),
        Comment                      => $comment,
    };
}

@out_records = sort {
       chrom_cmp($a->{Chrom}, $b->{Chrom})
    || $a->{Split_Start} <=> $b->{Split_Start}
    || $a->{Split_End}   <=> $b->{Split_End}
    || $a->{Cluster_ID} cmp $b->{Cluster_ID}
} @out_records;

write_output($opt{out}, \@out_records);

my $n_high = scalar(grep { $_->{Evidence_Level} eq 'High' } @out_records);
my $n_mod  = scalar(grep { $_->{Evidence_Level} eq 'Moderate' } @out_records);
my $n_low  = scalar(grep { $_->{Evidence_Level} eq 'Low' } @out_records);

log_msg("Output candidates written: " . scalar(@out_records));
log_msg("High=$n_high Moderate=$n_mod Low=$n_low");
log_msg("Merge evidence finished");

exit 0;

# ============================================================
# Readers
# ============================================================

sub read_split_records {
    my ($file, $sample_filter) = @_;

    my ($cols, $rows) = read_table($file);
    my %idx = index_cols($cols);

    my $sample_col = require_any_col($file, \%idx, qw/Sample SampleID/);
    my $id_col     = require_any_col($file, \%idx, qw/Cluster_ID ID/);
    my $chrom_col  = require_any_col($file, \%idx, qw/SV_Chrom Chrom CHROM/);
    my $start_col  = require_any_col($file, \%idx, qw/Cluster_Start Start START/);
    my $end_col    = require_any_col($file, \%idx, qw/Cluster_End End END/);

    my $gene_col   = find_col(\%idx, qw/Scan_Genes Gene Genes Target_Gene Target_Name/);
    my $status_col = find_col(\%idx, qw/Cluster_Status Status Candidate_Status/);

    my @out;

    for my $row (@$rows) {
        my $sample = normalize_value($row->{$sample_col});

        next if $sample ne $sample_filter;

        if (defined $status_col) {
            my $status = normalize_value($row->{$status_col});

            next if $status ne 'NA' && $status !~ /^Pass$/i;
        }

        my ($s, $e) = normalize_interval($row->{$start_col}, $row->{$end_col}, $file);

        my $gene = defined $gene_col ? normalize_value($row->{$gene_col}) : 'NA';

        push @out, {
            Type      => 'Split',
            SampleID  => $sample,
            Gene      => $gene,
            Chrom     => normalize_chrom_value($row->{$chrom_col}),
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

    my ($cols, $rows) = read_table($file);
    my %idx = index_cols($cols);

    my $sample_col     = find_col(\%idx, qw/SampleID Sample/);
    my $id_col         = find_col(\%idx, qw/Cluster_ID Del_Window_IDs Depth_File ID/);
    my $chrom_col      = require_any_col($file, \%idx, qw/Chrom CHROM SV_Chrom/);
    my $start_col      = require_any_col($file, \%idx, qw/Start START Cluster_Start/);
    my $end_col        = require_any_col($file, \%idx, qw/End END Cluster_End/);
    my $gene_col       = find_col(\%idx, qw/Gene Scan_Genes Genes Target_Gene Target_Name/);
    my $status_col     = find_col(\%idx, qw/Discordant_Status Candidate_Status Status Cluster_Status/);
    my $type_col       = find_col(\%idx, qw/Type TYPE SV_Type SVTYPE/);
    my $depth_file_col = find_col(\%idx, qw/Depth_File DepthFile File/);

    my @out;
    my $n = 0;

    for my $row (@$rows) {
        $n++;

        if (defined $sample_col) {
            my $sample = normalize_value($row->{$sample_col});
            next if $sample ne $sample_filter;
        }

        if ($type eq 'Depth' && defined $type_col) {
            my $svtype = normalize_value($row->{$type_col});
            next if $svtype ne 'NA' && $svtype !~ /DEL/i;
        }

        if ($type eq 'Discordant' && defined $status_col) {
            my $status = normalize_value($row->{$status_col});

            next if $status ne 'NA' && $status =~ /(discard|fail|low)/i;
        }

        my ($s, $e) = normalize_interval($row->{$start_col}, $row->{$end_col}, $file);

        my $gene = defined $gene_col ? normalize_value($row->{$gene_col}) : 'NA';

        if ($gene eq 'NA' && defined $depth_file_col) {
            $gene = guess_gene_from_depth_file($row->{$depth_file_col});
        }

        my $id = defined $id_col ? normalize_value($row->{$id_col}) : 'NA';
        $id = $type . '_' . $n if $id eq 'NA';

        push @out, {
            Type      => $type,
            SampleID  => $sample_filter,
            Gene      => $gene,
            Chrom     => normalize_chrom_value($row->{$chrom_col}),
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

    open my $fh, '<', $file
        or die "[ERROR] Cannot open $file: $!\n";

    my $header = <$fh>;

    die "[ERROR] Empty file: $file\n"
        unless defined $header;

    chomp $header;
    $header =~ s/\r$//;
    $header =~ s/^\s+|\s+$//g;

    my @cols = split_line($header);

    @cols = map {
        my $x = $_;
        $x =~ s/^\s+|\s+$//g;
        $x;
    } @cols;

    die "[ERROR] No header columns found in $file\n"
        unless @cols;

    for my $c (@cols) {
        die "[ERROR] Empty column name found in header of $file\n"
            unless defined $c && $c ne '';
    }

    my @rows;

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;

        my @v = split_line($line, scalar(@cols));

        my %row;

        for my $i (0 .. $#cols) {
            my $val = defined $v[$i] ? $v[$i] : 'NA';
            $val =~ s/^\s+|\s+$//g;
            $row{$cols[$i]} = $val ne '' ? $val : 'NA';
        }

        push @rows, \%row;
    }

    close $fh;

    return (\@cols, \@rows);
}

sub split_line {
    my ($line, $limit) = @_;

    $line = '' unless defined $line;

    if ($line =~ /\t/) {
        return split /\t/, $line, defined $limit ? $limit : 0;
    }

    $line =~ s/^\s+|\s+$//g;

    return split /\s+/, $line, defined $limit ? $limit : 0;
}

# ============================================================
# Evidence logic
# ============================================================

sub covered_bases_by_union {
    my ($records_ref, $query_start, $query_end) = @_;

    my @iv;

    for my $r (@$records_ref) {
        my $s = max_num($query_start, $r->{Start});
        my $e = min_num($query_end,   $r->{End});

        next if $s > $e;

        push @iv, [ $s, $e ];
    }

    return (0, 0) unless @iv;

    @iv = sort {
           $a->[0] <=> $b->[0]
        || $a->[1] <=> $b->[1]
    } @iv;

    my $covered = 0;
    my ($cs, $ce) = @{ shift @iv };

    for my $x (@iv) {
        my ($s, $e) = @$x;

        if ($s <= $ce + 1) {
            $ce = $e if $e > $ce;
        }
        else {
            $covered += $ce - $cs + 1;
            ($cs, $ce) = ($s, $e);
        }
    }

    $covered += $ce - $cs + 1;

    my $query_len = interval_len($query_start, $query_end);
    my $frac = $query_len > 0 ? $covered / $query_len : 0;

    return ($covered, $frac);
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

sub output_header {
    return qw/
        SampleID
        Gene
        Chrom
        Cluster_ID
        Split_Start
        Split_End
        Split_Size
        Split_Reads
        Split_Records
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
}

sub write_output {
    my ($file, $records_ref) = @_;

    my $dir = dirname($file);

    make_path($dir)
        if defined $dir && $dir ne '.' && !-d $dir;

    my @header = output_header();

    open my $out, '>', $file
        or die "[ERROR] Cannot write $file: $!\n";

    print $out join("\t", @header), "\n";

    for my $r (@$records_ref) {
        my @v = map {
            clean_tsv_value(
                defined $r->{$_} && $r->{$_} ne '' ? $r->{$_} : 'NA'
            )
        } @header;

        print $out join("\t", @v), "\n";
    }

    close $out;
}

sub make_comment {
    my ($level, $depth_frac, $discordant_frac, $threshold) = @_;

    $depth_frac      = sprintf('%.4f', $depth_frac);
    $discordant_frac = sprintf('%.4f', $discordant_frac);

    if ($level eq 'High') {
        return "Split candidate supported by depth and discordant evidence";
    }

    if ($level eq 'Moderate') {
        return "Split candidate supported by depth only; no discordant interval covers >= $threshold of split interval";
    }

    return "Split candidate lacks >= $threshold depth coverage; confidence forced to Low";
}

sub candidate_status {
    my ($level) = @_;

    return 'High_confidence_candidate'     if $level eq 'High';
    return 'Moderate_confidence_candidate' if $level eq 'Moderate';

    return 'Low_confidence_candidate';
}

sub ranges_for_records {
    my @records = @_;

    return 'NA' unless @records;

    return uniq_values(
        map {
            $_->{Chrom} . ':' . $_->{Start} . '-' . $_->{End}
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
        ':',
        $r->{Type},
        $r->{Source_ID},
        $r->{Chrom} . ':' . $r->{Start} . '-' . $r->{End}
    );
}

sub sum_meta_numeric {
    my ($records_ref, $field) = @_;

    my $sum = 0;
    my $has = 0;
    my %seen;

    for my $r (@$records_ref) {
        next if $seen{$r->{Source_ID}}++;

        my $v = meta_value($r, $field);

        next unless defined $v && $v =~ /^-?\d+(?:\.\d+)?$/;

        $sum += $v;
        $has = 1;
    }

    return $has ? $sum : 'NA';
}

sub meta_value {
    my ($r, $field) = @_;

    return 'NA'
        unless defined $r && exists $r->{Meta};

    return normalize_value($r->{Meta}{$field})
        if exists $r->{Meta}{$field};

    my $lc = lc $field;

    for my $k (keys %{ $r->{Meta} }) {
        return normalize_value($r->{Meta}{$k})
            if lc($k) eq $lc;
    }

    return 'NA';
}

sub clean_tsv_value {
    my ($v) = @_;

    $v = 'NA' unless defined $v;
    $v = 'NA' if $v eq '';

    $v =~ s/\r/ /g;
    $v =~ s/\n/ /g;
    $v =~ s/\t/ /g;

    return $v;
}

# ============================================================
# General helpers
# ============================================================

sub check_existing_file {
    my ($file, $desc) = @_;

    die "[ERROR] $desc not found: $file\n"
        unless defined $file && -e $file;

    die "[ERROR] $desc is not a regular file: $file\n"
        unless -f $file;

    die "[ERROR] $desc is empty and has no header: $file\n"
        unless -s $file;

    return 1;
}

sub index_by_chrom {
    my @records = @_;

    my %h;

    for my $r (@records) {
        my @keys = chrom_keys($r->{Chrom});

        for my $k (@keys) {
            push @{ $h{$k} }, $r;
        }
    }

    return %h;
}

sub records_for_chrom {
    my ($index_ref, $chrom) = @_;

    my @out;
    my %seen_ref;

    for my $k (chrom_keys($chrom)) {
        next unless exists $index_ref->{$k};

        for my $r (@{ $index_ref->{$k} }) {
            my $id = "$r";
            next if $seen_ref{$id}++;

            push @out, $r;
        }
    }

    return @out;
}

sub chrom_keys {
    my ($chrom) = @_;

    $chrom = normalize_chrom_value($chrom);

    return ('NA') if $chrom eq 'NA';

    my @keys;
    my %seen;

    my $add = sub {
        my ($x) = @_;
        return unless defined $x && $x ne '';
        return if $seen{$x}++;
        push @keys, $x;
    };

    $add->($chrom);

    if ($chrom =~ /^chr(.+)$/i) {
        my $no_chr = $1;
        $add->($no_chr);

        if (uc($no_chr) eq 'M') {
            $add->('MT');
            $add->('chrMT');
        }
        elsif (uc($no_chr) eq 'MT') {
            $add->('M');
            $add->('chrM');
        }
    }
    else {
        $add->("chr$chrom");

        if (uc($chrom) eq 'M') {
            $add->('MT');
            $add->('chrMT');
        }
        elsif (uc($chrom) eq 'MT') {
            $add->('M');
            $add->('chrM');
        }
    }

    return @keys;
}

sub normalize_chrom_value {
    my ($v) = @_;

    $v = normalize_value($v);

    return 'NA' if $v eq 'NA';

    $v =~ s/^\s+|\s+$//g;

    return $v;
}

sub index_cols {
    my ($cols_ref) = @_;

    my %idx;

    for my $c (@$cols_ref) {
        next unless defined $c;

        my $raw = $c;
        my $lc  = lc $c;

        die "[ERROR] Duplicate column name found: $raw\n"
            if exists $idx{$raw};

        $idx{$raw} = $raw;
        $idx{$lc}  = $raw;
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
      . join(',', @candidates)
      . "\n"
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

    return 'NA' unless defined $v;

    $v =~ s/^\s+//;
    $v =~ s/\s+$//;

    return $v eq '' ? 'NA' : $v;
}

sub guess_gene_from_depth_file {
    my ($v) = @_;

    $v = normalize_value($v);

    return 'NA' if $v eq 'NA';

    my $base = $v;

    $base =~ s{.*/}{};
    $base =~ s/\.depth\.tsv$//;

    if ($base =~ /_([A-Za-z0-9.-]+)_N[MR]_\d+/) {
        return $1;
    }

    if ($base =~ /^HCM_([A-Za-z0-9.-]+)_DEL/) {
        return $1;
    }

    return 'NA';
}

sub uniq_values {
    my @v = @_;

    my %seen;
    my @out;

    for my $x (@v) {
        $x = normalize_value($x);

        next if $x eq 'NA';
        next if $seen{$x}++;

        push @out, $x;
    }

    return @out ? join(',', @out) : 'NA';
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

    return 23 if uc($c) eq 'X';
    return 24 if uc($c) eq 'Y';
    return 25 if uc($c) =~ /^(M|MT)$/;

    return 1000;
}

sub min_num {
    return $_[0] < $_[1] ? $_[0] : $_[1];
}

sub max_num {
    return $_[0] > $_[1] ? $_[0] : $_[1];
}

sub log_msg {
    my ($msg) = @_;

    print STDERR "[INFO] $msg\n";
}

sub usage {
    return <<'USAGE';
Usage:
  perl bin/merge_evidence.pl \
    --sample SAMPLE_ID \
    --depth SAMPLE.depth_candidates.cusum.all.tsv \
    --split SAMPLE.split_reads.clusters.tsv \
    --discordant SAMPLE.discordant_reads.tsv \
    --out SAMPLE.merged_candidates.tsv

Optional:
  --min-overlap FLOAT
      Minimum fraction of split interval covered by supporting evidence.
      Default: 0.90

Merging logic:
  1. Use split_reads.clusters.tsv as the core candidate list.
  2. For each passed split cluster, calculate how many split bases are
     covered by overlapping CUSUM depth intervals.
  3. Depth intervals are merged by union before coverage calculation.
  4. For discordant evidence, require at least one discordant-read cluster
     to cover >= min-overlap of the split interval.
  5. Confidence:
       High     = Split + Depth + Discordant
       Moderate = Split + Depth only
       Low      = Split without Depth

Empty-result behavior:
  If no passed split-read cluster exists, write a header-only output file
  and exit normally.

Coordinates:
  1-based closed intervals.
USAGE
}

