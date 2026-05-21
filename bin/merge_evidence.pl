use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename qw(dirname);

# merge_evidence.pl
#
# Purpose:
#   Merge exon-level deletion evidence from three independent modules:
#
#     1) Depth evidence
#        Input: depth_candidates.tsv
#        Meaning: continuous low-depth windows suggesting deletion.
#        Limitation: boundary is affected by window size and is usually rough.
#
#     2) Split-read evidence
#        Input: split_reads.clusters.tsv
#        Meaning: split reads / SA-tag reads supporting breakpoint-spanning events.
#        Strength: usually closest to the real breakpoint.
#
#     3) Discordant-read evidence
#        Input: discordant_reads.tsv
#        Meaning: read pairs with abnormal insert size / orientation.
#        Strength: supports deletion-like SV, but breakpoint is usually less precise
#                  than split-read evidence.
#
# Main merging logic:
#
#   For each Gene + Chrom group:
#
#     Step 1. Convert all evidence intervals into position-level support events.
#             This is implemented by interval event scanning, not by expanding
#             every single base.
#
#     Step 2. For each continuous genomic segment, count whether it is supported by:
#
#               Depth
#               Split
#               Discordant
#
#             Evidence_Count = number of supporting evidence types.
#
#     Step 3. The main output keeps only regions with:
#
#               Evidence_Count == 3
#
#             That means the final main candidates are supported by:
#
#               Depth + Split + Discordant
#
#     Step 4. An auxiliary file <out>.all.tsv is also written.
#             It keeps all regions with Evidence_Count >= 1 for debugging.
#
#
# Coordinate definitions:
#
#   This script outputs three coordinate systems:
#
#   1) Merged_Start / Merged_End
#
#      Definition:
#        The outer boundary of all evidence intervals overlapping this segment.
#
#      Calculation:
#        Merged_Start = minimum start among overlapping evidence records
#        Merged_End   = maximum end among overlapping evidence records
#
#      Meaning:
#        The widest possible event range supported by any evidence.
#        This is useful for checking the complete evidence span, but it may be
#        wider than the true deletion boundary.
#
#
#   2) Core_Start / Core_End
#
#      Definition:
#        The position-level segment currently supported by the listed evidence types.
#
#      For three-evidence candidates, this is the intersection of:
#
#        Depth interval
#        Split interval
#        Discordant interval
#
#      Meaning:
#        The most conservative region jointly supported by all available evidence.
#        It is useful as a high-confidence core region, but it may be shorter than
#        the real deletion, because depth windows are often shorter or coarser.
#
#
#   3) Best_Start / Best_End
#
#      Definition:
#        The recommended candidate boundary for downstream annotation and reporting.
#
#      Priority:
#
#        Split > Discordant > Depth
#
#      Rationale:
#        Split-read clusters are usually closest to the real breakpoint.
#        Discordant-read clusters are the second choice.
#        Depth intervals are used only when no breakpoint-level evidence exists.
#
#      For example:
#
#        Depth      chr18:34232176-34241175
#        Split      chr18:34232240-34241309
#        Discordant chr18:34232175-34241310
#
#        Merged     chr18:34232175-34241310
#        Core       chr18:34232240-34241175
#        Best       chr18:34232240-34241309   derived from Split
#
#
# Input requirements:
#
#   --depth:
#     depth_candidates.tsv
#     Required columns:
#       SampleID Gene Chrom Start End
#
#   --split:
#     split_reads.clusters.tsv only
#     Required columns:
#       Sample Cluster_ID SV_Chrom Cluster_Start Cluster_End
#
#   --discordant:
#     discordant_reads.tsv
#     Required columns:
#       SampleID Cluster_ID Chrom Start End
#
# Output:
#
#   --out:
#     Main result. Only three-evidence candidates.
#
#   <out>.all.tsv:
#     Auxiliary result. All evidence-supported segments.
#
# Coordinates:
#   All coordinates are treated as 1-based closed intervals.

my %opt;

GetOptions(
    'config=s'     => \$opt{config},
    'sample=s'     => \$opt{sample},
    'depth=s'      => \$opt{depth},
    'split=s'      => \$opt{split},
    'discordant=s' => \$opt{discordant},
    'out=s'        => \$opt{out},
    'help'         => \$opt{help},
) or die usage();

if ($opt{help}) {
    print usage();
    exit 0;
}

for my $k (qw/config sample depth split discordant out/) {
    die "[ERROR] Missing required option --$k\n" . usage()
        unless defined $opt{$k} && $opt{$k} ne '';
}

die "[ERROR] Config file not found: $opt{config}\n"         unless -s $opt{config};
die "[ERROR] Depth file not found: $opt{depth}\n"           unless -s $opt{depth};
die "[ERROR] Split cluster file not found: $opt{split}\n"   unless -s $opt{split};
die "[ERROR] Discordant file not found: $opt{discordant}\n" unless -s $opt{discordant};

if ($opt{split} =~ /split_reads\.tsv$/ && $opt{split} !~ /clusters/) {
    die "[ERROR] Wrong split input file: $opt{split}\n"
      . "        merge_evidence.pl requires split_reads.clusters.tsv, not split_reads.tsv\n";
}

my $all_out = $opt{out} . '.all.tsv';

log_msg("Merge evidence started");
log_msg("Config          : $opt{config}");
log_msg("Sample          : $opt{sample}");
log_msg("Depth file      : $opt{depth}");
log_msg("Split file      : $opt{split}");
log_msg("Discordant file : $opt{discordant}");
log_msg("Main output     : $opt{out}");
log_msg("All output      : $all_out");

my @records;

push @records, read_depth_records($opt{depth}, $opt{sample});
push @records, read_split_cluster_records($opt{split}, $opt{sample});
push @records, read_discordant_records($opt{discordant}, $opt{sample});

die "[ERROR] No evidence records loaded from input files\n" unless @records;

log_msg("Total evidence records loaded: " . scalar(@records));

my %groups;

for my $r (@records) {
    my $gene  = normalize_value($r->{Gene});
    my $chrom = normalize_value($r->{Chrom});

    next if $gene eq 'NA' || $chrom eq 'NA';

    my $key = join("\t", $gene, $chrom);
    push @{ $groups{$key} }, $r;
}

my @all_candidates;

for my $key (sort keys %groups) {
    my ($gene, $chrom) = split /\t/, $key;

    my @segments = build_position_support_segments(
        $gene,
        $chrom,
        $groups{$key}
    );

    push @all_candidates, @segments;
}

@all_candidates = sort {
       chrom_cmp($a->{Chrom}, $b->{Chrom})
    || $a->{Core_Start} <=> $b->{Core_Start}
    || $a->{Core_End}   <=> $b->{Core_End}
    || $a->{Gene} cmp $b->{Gene}
} @all_candidates;

my @main_candidates = grep {
    $_->{Evidence_Count} == 3
} @all_candidates;

write_candidates($all_out, \@all_candidates);
write_candidates($opt{out}, \@main_candidates);

log_msg("All evidence regions written       : " . scalar(@all_candidates));
log_msg("Three-evidence candidates written  : " . scalar(@main_candidates));
log_msg("Merge evidence finished");

exit 0;


sub usage {
    return <<'USAGE';

Usage:
  perl merge_evidence.pl \
    --config conf/hcm_exondel.example.conf \
    --sample 25B09089386 \
    --depth 01.depth/25B09089386.depth_candidates.tsv \
    --split 02.split_reads/25B09089386.split_reads.clusters.tsv \
    --discordant 03.discordant_reads/25B09089386.discordant_reads.tsv \
    --out 04.candidates/25B09089386.merged_candidates.tsv

Required:
  --config       Config file. Required for compatibility with main pipeline.
  --sample       Sample ID.
  --depth        Depth candidate file.
  --split        Split-read cluster file. Only split_reads.clusters.tsv is supported.
  --discordant   Discordant-read candidate file.
  --out          Main output file.

Output:
  --out          Main output. Only Depth + Split + Discordant regions.
  <out>.all.tsv  Auxiliary output. All evidence-supported regions.

Coordinate meaning:
  Merged_*       Outer boundary of overlapping evidence intervals.
  Core_*         Conservative position-level evidence-overlap region.
  Best_*         Recommended boundary for annotation/reporting.
                 Priority: Split > Discordant > Depth.

Notes:
  Coordinates are treated as 1-based closed intervals.

USAGE
}


sub log_msg {
    my ($msg) = @_;
    print STDERR "[INFO] $msg\n";
}


sub read_table {
    my ($file) = @_;

    die "[ERROR] File not found or empty: $file\n" unless -s $file;

    open my $fh, '<', $file
        or die "[ERROR] Cannot open file $file: $!\n";

    my $header = <$fh>;
    die "[ERROR] Empty file: $file\n" unless defined $header;

    chomp $header;
    $header =~ s/\r$//;

    my @cols = split_line($header);
    my @rows;

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;

        my @v = split_line($line, scalar(@cols));

        my %row;
        for my $i (0 .. $#cols) {
            $row{$cols[$i]} =
                defined $v[$i] && $v[$i] ne ''
                ? $v[$i]
                : 'NA';
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
    else {
        return split /\s+/, $line, defined $limit ? $limit : 0;
    }
}


sub require_columns {
    my ($file, $cols_ref, @required) = @_;

    my %has = map { $_ => 1 } @$cols_ref;

    for my $c (@required) {
        die "[ERROR] Required column '$c' not found in $file\n"
            unless $has{$c};
    }
}


sub read_depth_records {
    my ($file, $sample_filter) = @_;

    my ($cols, $rows) = read_table($file);

    require_columns(
        $file,
        $cols,
        qw/SampleID Gene Chrom Start End/
    );

    my @out;

    for my $row (@$rows) {
        next if normalize_value($row->{SampleID}) ne $sample_filter;

        my ($s, $e) = normalize_interval(
            $row->{Start},
            $row->{End},
            $file
        );

        push @out, {
            Type           => 'Depth',
            SampleID       => normalize_value($row->{SampleID}),
            Gene           => normalize_value($row->{Gene}),
            Classification => normalize_value($row->{Classification}),
            Transcript     => normalize_value($row->{Transcript}),
            Chrom          => normalize_value($row->{Chrom}),
            Start          => $s,
            End            => $e,
            Source_ID      => normalize_value($row->{Del_Window_IDs}),
            Meta           => {%$row},
        };
    }

    log_msg("Depth records loaded: " . scalar(@out));

    return @out;
}


sub read_split_cluster_records {
    my ($file, $sample_filter) = @_;

    my ($cols, $rows) = read_table($file);

    require_columns(
        $file,
        $cols,
        qw/Sample Cluster_ID SV_Chrom Cluster_Start Cluster_End/
    );

    my @out;

    for my $row (@$rows) {
        next if normalize_value($row->{Sample}) ne $sample_filter;

        if (exists $row->{Cluster_Status}) {
            next if normalize_value($row->{Cluster_Status}) ne 'Pass';
        }

        my ($s, $e) = normalize_interval(
            $row->{Cluster_Start},
            $row->{Cluster_End},
            $file
        );

        my @genes = split_genes(first_non_na(
            $row->{Scan_Genes},
            $row->{Gene},
            $row->{Genes}
        ));

        @genes = ('NA') unless @genes;

        for my $gene (@genes) {
            push @out, {
                Type           => 'Split',
                SampleID       => normalize_value($row->{Sample}),
                Gene           => $gene,
                Classification => 'NA',
                Transcript     => 'NA',
                Chrom          => normalize_value($row->{SV_Chrom}),
                Start          => $s,
                End            => $e,
                Source_ID      => normalize_value($row->{Cluster_ID}),
                Meta           => {%$row},
            };
        }
    }

    log_msg("Split cluster records loaded: " . scalar(@out));

    return @out;
}


sub read_discordant_records {
    my ($file, $sample_filter) = @_;

    my ($cols, $rows) = read_table($file);

    require_columns(
        $file,
        $cols,
        qw/SampleID Cluster_ID Chrom Start End/
    );

    my @out;

    for my $row (@$rows) {
        next if normalize_value($row->{SampleID}) ne $sample_filter;

        my ($s, $e) = normalize_interval(
            $row->{Start},
            $row->{End},
            $file
        );

        my @genes = split_genes(first_non_na(
            $row->{Gene},
            $row->{Target_Gene},
            $row->{Target_Name}
        ));

        @genes = ('NA') unless @genes;

        for my $gene (@genes) {
            push @out, {
                Type           => 'Discordant',
                SampleID       => normalize_value($row->{SampleID}),
                Gene           => $gene,
                Classification => 'NA',
                Transcript     => normalize_value($row->{Transcript}),
                Chrom          => normalize_value($row->{Chrom}),
                Start          => $s,
                End            => $e,
                Source_ID      => normalize_value($row->{Cluster_ID}),
                Meta           => {%$row},
            };
        }
    }

    log_msg("Discordant records loaded: " . scalar(@out));

    return @out;
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


sub first_non_na {
    for my $v (@_) {
        $v = normalize_value($v);
        return $v if $v ne 'NA';
    }

    return 'NA';
}


sub split_genes {
    my ($v) = @_;

    $v = normalize_value($v);

    return () if $v eq 'NA' || $v eq '.' || $v eq '-';

    if ($v =~ /^([A-Za-z0-9_.-]+)\|/) {
        $v = $1;
    }

    my @genes = split /[,;]+/, $v;

    @genes = grep {
        $_ ne '' && $_ ne 'NA'
    } map {
        normalize_value($_)
    } @genes;

    my %seen;

    return grep {
        !$seen{$_}++
    } @genes;
}


sub build_position_support_segments {
    my ($gene, $chrom, $records) = @_;

    my %events;

    for my $r (@$records) {
        my $type = $r->{Type};

        $events{$r->{Start}}{$type}++;
        $events{$r->{End} + 1}{$type}--;
    }

    my @pos = sort { $a <=> $b } keys %events;

    return () if @pos < 2;

    my %active = (
        Depth      => 0,
        Split      => 0,
        Discordant => 0,
    );

    my @segments;

    for my $i (0 .. $#pos - 1) {
        my $p = $pos[$i];

        for my $type (keys %{ $events{$p} }) {
            $active{$type} += $events{$p}{$type};
        }

        my $seg_start = $p;
        my $seg_end   = $pos[$i + 1] - 1;

        next if $seg_start > $seg_end;

        my @types = grep {
            $active{$_} > 0
        } qw/Depth Split Discordant/;

        next unless @types;

        my @overlap_records = grep {
            $_->{Start} <= $seg_end && $_->{End} >= $seg_start
        } @$records;

        push @segments, make_candidate(
            $gene,
            $chrom,
            $seg_start,
            $seg_end,
            \@types,
            \@overlap_records
        );
    }

    @segments = merge_adjacent_segments(@segments);

    return @segments;
}


sub make_candidate {
    my ($gene, $chrom, $core_start, $core_end, $types_ref, $records_ref) = @_;

    my %type_has = map { $_ => 1 } @$types_ref;

    my @ordered_types = grep {
        $type_has{$_}
    } qw/Depth Split Discordant/;

    my $evidence_count = scalar(@ordered_types);

    my ($merged_start, $merged_end) = ($core_start, $core_end);

    for my $r (@$records_ref) {
        $merged_start = $r->{Start} if $r->{Start} < $merged_start;
        $merged_end   = $r->{End}   if $r->{End}   > $merged_end;
    }

    my @depth_records      = grep { $_->{Type} eq 'Depth' } @$records_ref;
    my @split_records      = grep { $_->{Type} eq 'Split' } @$records_ref;
    my @discordant_records = grep { $_->{Type} eq 'Discordant' } @$records_ref;

    my ($best_start, $best_end, $best_evidence) = select_best_interval(
        \@depth_records,
        \@split_records,
        \@discordant_records
    );

    my $best_size = 'NA';
    if ($best_start ne 'NA' && $best_end ne 'NA') {
        $best_size = $best_end - $best_start + 1;
    }

    my $sample = first_non_na(
        (map { $_->{SampleID} } @$records_ref)
    );

    my $classification = first_non_na(
        (map { $_->{Classification} } @depth_records)
    );

    my $transcript = first_non_na(
        (map { $_->{Transcript} } @depth_records),
        (map { $_->{Transcript} } @discordant_records)
    );

    my $evidence_types = join(',', @ordered_types);
    my $level = evidence_level($evidence_count);

    return {
        SampleID                => $sample,
        Gene                    => $gene,
        Classification          => $classification,
        Transcript              => $transcript,
        Chrom                   => $chrom,

        Merged_Start            => $merged_start,
        Merged_End              => $merged_end,
        Merged_Size             => $merged_end - $merged_start + 1,

        Core_Start              => $core_start,
        Core_End                => $core_end,
        Core_Size               => $core_end - $core_start + 1,

        Best_Start              => $best_start,
        Best_End                => $best_end,
        Best_Size               => $best_size,
        Best_Evidence           => $best_evidence,

        Evidence_Count          => $evidence_count,
        Evidence_Level          => $level,
        Evidence_Types          => $evidence_types,

        Depth_Support           => $type_has{Depth}      ? 'Yes' : 'No',
        Split_Read_Support      => $type_has{Split}      ? 'Yes' : 'No',
        Discordant_Read_Support => $type_has{Discordant} ? 'Yes' : 'No',

        Depth_Range             => ranges_for_type(\@depth_records),
        Split_Range             => ranges_for_type(\@split_records),
        Discordant_Range        => ranges_for_type(\@discordant_records),

        Window_Count            => sum_unique_field(\@depth_records, 'Source_ID', 'Window_Count'),
        Gene_Mean_Depth         => uniq_field(\@depth_records, 'Gene_Mean_Depth'),
        Candidate_Mean_Depth    => uniq_field(\@depth_records, 'Candidate_Mean_Depth'),
        Mean_Depth_Ratio        => uniq_field(\@depth_records, 'Mean_Depth_Ratio'),
        Min_Depth_Ratio         => uniq_field(\@depth_records, 'Min_Depth_Ratio'),
        Max_Depth_Ratio         => uniq_field(\@depth_records, 'Max_Depth_Ratio'),

        Split_Cluster_IDs       => uniq_values((map { $_->{Source_ID} } @split_records)),
        Split_Reads             => sum_unique_field(\@split_records, 'Source_ID', 'Support_Reads'),
        Split_Records           => sum_unique_field(\@split_records, 'Source_ID', 'Support_Records'),

        Discordant_Cluster_IDs  => uniq_values((map { $_->{Source_ID} } @discordant_records)),
        Discordant_Reads        => sum_unique_field(\@discordant_records, 'Source_ID', 'Discordant_Reads'),
        Median_Insert_Size      => uniq_field(\@discordant_records, 'Median_Insert_Size'),

        All_Genes               => uniq_values((map { $_->{Gene} } @$records_ref)),
        All_Transcripts         => uniq_values((map { $_->{Transcript} } @$records_ref)),
        Source_Records          => source_records($records_ref),

        Candidate_Status        => $evidence_count == 3 ? 'Three_evidence_candidate'
                                  : $evidence_count == 2 ? 'Two_evidence_candidate'
                                  : 'Single_evidence_candidate',

        Comment                 => "Core region supported by $evidence_types evidence; Best interval selected by priority Split > Discordant > Depth",
    };
}


sub select_best_interval {
    my ($depth_records, $split_records, $discordant_records) = @_;

    if (@$split_records) {
        return representative_interval($split_records, 'Split');
    }

    if (@$discordant_records) {
        return representative_interval($discordant_records, 'Discordant');
    }

    if (@$depth_records) {
        return representative_interval($depth_records, 'Depth');
    }

    return ('NA', 'NA', 'NA');
}


sub representative_interval {
    my ($records_ref, $evidence_name) = @_;

    return ('NA', 'NA', 'NA') unless @$records_ref;

    my $start = $records_ref->[0]{Start};
    my $end   = $records_ref->[0]{End};

    for my $r (@$records_ref) {
        $start = $r->{Start} if $r->{Start} < $start;
        $end   = $r->{End}   if $r->{End}   > $end;
    }

    return ($start, $end, $evidence_name);
}


sub merge_adjacent_segments {
    my @segments = @_;

    return @segments unless @segments;

    @segments = sort {
           $a->{Core_Start} <=> $b->{Core_Start}
        || $a->{Core_End}   <=> $b->{Core_End}
    } @segments;

    my @merged;

    for my $s (@segments) {
        if (
            @merged
            && $merged[-1]{Chrom} eq $s->{Chrom}
            && $merged[-1]{Gene} eq $s->{Gene}
            && $merged[-1]{Evidence_Types} eq $s->{Evidence_Types}
            && $merged[-1]{Best_Evidence} eq $s->{Best_Evidence}
            && $merged[-1]{Core_End} + 1 == $s->{Core_Start}
        ) {
            my $prev = $merged[-1];

            $prev->{Core_End}  = $s->{Core_End};
            $prev->{Core_Size} = $prev->{Core_End} - $prev->{Core_Start} + 1;

            $prev->{Merged_Start} = min_num($prev->{Merged_Start}, $s->{Merged_Start});
            $prev->{Merged_End}   = max_num($prev->{Merged_End},   $s->{Merged_End});
            $prev->{Merged_Size}  = $prev->{Merged_End} - $prev->{Merged_Start} + 1;

            if ($prev->{Best_Start} ne 'NA' && $s->{Best_Start} ne 'NA') {
                $prev->{Best_Start} = min_num($prev->{Best_Start}, $s->{Best_Start});
                $prev->{Best_End}   = max_num($prev->{Best_End},   $s->{Best_End});
                $prev->{Best_Size}  = $prev->{Best_End} - $prev->{Best_Start} + 1;
            }

            for my $field (
                qw/
                    Depth_Range
                    Split_Range
                    Discordant_Range
                    Split_Cluster_IDs
                    Discordant_Cluster_IDs
                    All_Genes
                    All_Transcripts
                    Source_Records
                /
            ) {
                $prev->{$field} = merge_csv_values(
                    $prev->{$field},
                    $s->{$field}
                );
            }
        }
        else {
            push @merged, {%$s};
        }
    }

    return @merged;
}


sub ranges_for_type {
    my ($records_ref) = @_;

    return 'NA' unless @$records_ref;

    return uniq_values(
        (map {
            $_->{Chrom} . ':' . $_->{Start} . '-' . $_->{End}
        } @$records_ref)
    );
}


sub source_records {
    my ($records_ref) = @_;

    return uniq_values(
        (map {
            $_->{Type}
            . ':'
            . $_->{Source_ID}
            . ':'
            . $_->{Chrom}
            . ':'
            . $_->{Start}
            . '-'
            . $_->{End}
        } @$records_ref)
    );
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


sub uniq_field {
    my ($records_ref, $field) = @_;

    return uniq_values(
        (map {
            $_->{Meta}{$field}
        } @$records_ref)
    );
}


sub sum_unique_field {
    my ($records_ref, $id_field, $value_field) = @_;

    my %seen;
    my $sum = 0;
    my $has = 0;

    for my $r (@$records_ref) {
        my $id = normalize_value($r->{$id_field});
        $id = normalize_value($r->{Source_ID}) if $id eq 'NA';

        next if $seen{$id}++;

        my $v = $r->{Meta}{$value_field};

        next unless defined $v && $v =~ /^-?\d+(?:\.\d+)?$/;

        $sum += $v;
        $has = 1;
    }

    return $has ? $sum : 'NA';
}


sub merge_csv_values {
    my ($a, $b) = @_;

    my @v;

    for my $x ($a, $b) {
        next unless defined $x;
        next if $x eq 'NA';

        push @v, split /,/, $x;
    }

    return uniq_values(@v);
}


sub evidence_level {
    my ($n) = @_;

    return $n >= 3 ? 'High'
         : $n == 2 ? 'Moderate'
         : 'Low';
}


sub write_candidates {
    my ($file, $candidates_ref) = @_;

    my $dir = dirname($file);

    if (defined $dir && $dir ne '.' && !-d $dir) {
        system('mkdir', '-p', $dir) == 0
            or die "[ERROR] Cannot create output directory: $dir\n";
    }

    my @header = qw/
        SampleID
        Gene
        Classification
        Transcript
        Chrom

        Merged_Start
        Merged_End
        Merged_Size

        Core_Start
        Core_End
        Core_Size

        Best_Start
        Best_End
        Best_Size
        Best_Evidence

        Evidence_Count
        Evidence_Level
        Evidence_Types

        Depth_Support
        Split_Read_Support
        Discordant_Read_Support

        Depth_Range
        Split_Range
        Discordant_Range

        Window_Count
        Gene_Mean_Depth
        Candidate_Mean_Depth
        Mean_Depth_Ratio
        Min_Depth_Ratio
        Max_Depth_Ratio

        Split_Cluster_IDs
        Split_Reads
        Split_Records

        Discordant_Cluster_IDs
        Discordant_Reads
        Median_Insert_Size

        All_Genes
        All_Transcripts
        Source_Records

        Candidate_Status
        Comment
    /;

    open my $out, '>', $file
        or die "[ERROR] Cannot write output file $file: $!\n";

    print $out join("\t", @header), "\n";

    for my $c (@$candidates_ref) {
        my @v = map {
            defined $c->{$_} && $c->{$_} ne '' ? $c->{$_} : 'NA'
        } @header;

        print $out join("\t", @v), "\n";
    }

    close $out;
}


sub chrom_cmp {
    my ($a, $b) = @_;

    my $aa = chrom_rank($a);
    my $bb = chrom_rank($b);

    return $aa <=> $bb || $a cmp $b;
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


