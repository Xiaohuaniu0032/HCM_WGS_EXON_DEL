#!/usr/bin/env perl
# -*- coding: utf-8 -*-

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path);

# ============================================================
# run_gene_mean_depth.pl
#
# Function:
#   Window-based depth-ratio analysis for exon-level deletion
#   detection in HCM-related genes.
#
# Main design:
#   1. Gene regions are read from REFSEQ_MANE_SELECT_GENE_TXT.
#   2. By default, only genes in HCM_CORE_GENE_LIST are analyzed.
#   3. Coordinates are 1-based closed intervals.
#   4. For each gene, samtools depth is called only once.
#   5. Gene mean depth and all window mean depths are calculated
#      from the per-gene base-level depth file.
#
# Required input:
#   --config
#   --bam
#   --sample
#   --out
#
# Required config:
#   SAMTOOLS
#   REFSEQ_MANE_SELECT_GENE_TXT
#   HCM_CORE_GENE_LIST
#
# Output:
#   1. *.gene_mean_depth.tsv
#   2. *.window_depth.tsv
#   3. *.all_window_ratio.tsv
#   4. *.depth_candidates.tsv
# ============================================================

my $config;
my $bam;
my $sample;
my $out;
my $threads = 4;
my $help = 0;

GetOptions(
    "config=s"  => \$config,
    "bam=s"     => \$bam,
    "sample=s"  => \$sample,
    "out=s"     => \$out,
    "threads=i" => \$threads,
    "help"      => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless $config && $bam && $sample && $out;

$config = abs_path($config);
$bam    = abs_path($bam);

die "[ERROR] Config file not found: $config\n" unless $config && -s $config;
die "[ERROR] BAM file not found: $bam\n" unless $bam && -s $bam;

my %CONF = read_config($config);

# ------------------------------------------------------------
# Config parameters
# ------------------------------------------------------------

my $samtools = get_conf(\%CONF, "SAMTOOLS", "samtools");

my $gene_txt = get_conf(\%CONF, "REFSEQ_MANE_SELECT_GENE_TXT", "");
$gene_txt = resolve_path_by_config($gene_txt, $config);

die "[ERROR] REFSEQ_MANE_SELECT_GENE_TXT file not found: $gene_txt\n"
    unless $gene_txt && -s $gene_txt;

my $analyze_core_only = get_conf(\%CONF, "ANALYZE_CORE_GENES_ONLY", 1);

my %core_gene;
my %gene_classification;
my $core_gene_list = get_conf(\%CONF, "HCM_CORE_GENE_LIST", "");

if ($analyze_core_only) {
    die "[ERROR] HCM_CORE_GENE_LIST is required when ANALYZE_CORE_GENES_ONLY=1\n"
        unless defined $core_gene_list && $core_gene_list ne "";

    $core_gene_list = resolve_path_by_config($core_gene_list, $config);

    die "[ERROR] HCM core gene list file not found: $core_gene_list\n"
        unless $core_gene_list && -s $core_gene_list;

    (%core_gene, %gene_classification) = read_core_gene_list($core_gene_list);
}

my $window_size = get_conf(\%CONF, "WINDOW_SIZE", 1000);
my $window_step = get_conf(\%CONF, "WINDOW_STEP", 500);
my $min_window_size = get_conf(\%CONF, "MIN_WINDOW_SIZE", 100);

my $min_gene_mean_depth   = get_conf(\%CONF, "MIN_GENE_MEAN_DEPTH", 20);
my $min_window_mean_depth = get_conf(\%CONF, "MIN_WINDOW_MEAN_DEPTH", 5);

my $min_depth_mapq = get_conf(\%CONF, "MIN_DEPTH_MAPQ", get_conf(\%CONF, "MIN_MAPQ", 20));
my $exclude_dup    = get_conf(\%CONF, "EXCLUDE_DUPLICATES_FOR_DEPTH", 1);

my $depth_ratio_cutoff = get_conf(\%CONF, "DEPTH_RATIO_CUTOFF", 0.65);
my $het_del_ratio_low  = get_conf(\%CONF, "HET_DEL_RATIO_LOW", 0.35);
my $het_del_ratio_high = get_conf(\%CONF, "HET_DEL_RATIO_HIGH", 0.70);

my $min_candidate_windows    = get_conf(\%CONF, "MIN_CANDIDATE_WINDOWS", 1);
my $max_merge_gap_windows    = get_conf(\%CONF, "MAX_MERGE_GAP_WINDOWS", 0);
my $min_depth_candidate_size = get_conf(\%CONF, "MIN_DEPTH_CANDIDATE_SIZE", 100);
my $max_depth_candidate_size = get_conf(\%CONF, "MAX_DEPTH_CANDIDATE_SIZE", 200000);

my $keep_depth_file = get_conf(\%CONF, "KEEP_GENE_DEPTH_FILE", 0);

check_positive_integer("WINDOW_SIZE", $window_size);
check_positive_integer("WINDOW_STEP", $window_step);
check_positive_integer("MIN_WINDOW_SIZE", $min_window_size);
check_non_negative_integer("MAX_MERGE_GAP_WINDOWS", $max_merge_gap_windows);

# ------------------------------------------------------------
# Output files
# ------------------------------------------------------------

my $outdir = dirname($out);
make_path($outdir) unless -d $outdir;

my $prefix = $out;
$prefix =~ s/\.tsv$//;

my $depth_dir = $prefix . ".gene_depth_files";
make_path($depth_dir) unless -d $depth_dir;

my $gene_depth_out   = $prefix . ".gene_mean_depth.tsv";
my $window_depth_out = $prefix . ".window_depth.tsv";
my $all_ratio_out    = $prefix . ".all_window_ratio.tsv";

# ============================================================
# Step 1. Read gene regions
# ============================================================

my @genes = read_gene_txt(
    file              => $gene_txt,
    analyze_core_only => $analyze_core_only,
    core_gene_ref     => \%core_gene,
    classification_ref => \%gene_classification,
);

die "[ERROR] No valid gene regions found in: $gene_txt\n" unless @genes;

print "[INFO] Sample: $sample\n";
print "[INFO] Gene coordinate file: $gene_txt\n";

if ($analyze_core_only) {
    print "[INFO] Core gene list: $core_gene_list\n";
    print "[INFO] Core genes loaded: " . scalar(keys %core_gene) . "\n";
    print "[INFO] Gene regions retained for analysis: " . scalar(@genes) . "\n";
} else {
    print "[INFO] ANALYZE_CORE_GENES_ONLY=0, all genes in gene TXT will be analyzed\n";
    print "[INFO] Gene regions retained for analysis: " . scalar(@genes) . "\n";
}

# ============================================================
# Step 2. Generate windows
# ============================================================

my @windows;
my %gene_windows;

foreach my $gene (@genes) {
    my $gkey = gene_key($gene);

    my @gene_windows = generate_windows(
        gene            => $gene,
        sample          => $sample,
        window_size     => $window_size,
        window_step     => $window_step,
        min_window_size => $min_window_size,
    );

    push @windows, @gene_windows;
    $gene_windows{$gkey} = \@gene_windows;
}

die "[ERROR] No windows generated. Please check gene TXT and window parameters\n"
    unless @windows;

# ============================================================
# Step 3. Generate per-gene depth file and calculate depths
# ============================================================

my %gene_mean_depth;
my %gene_total_length;
my %window_mean_depth;

foreach my $gene (@genes) {
    my $gkey = gene_key($gene);

    my $depth_file = make_gene_depth_filename(
        depth_dir => $depth_dir,
        sample    => $sample,
        gene      => $gene,
    );

    run_samtools_depth_for_gene(
        samtools    => $samtools,
        bam         => $bam,
        gene        => $gene,
        depth_file  => $depth_file,
        min_mapq    => $min_depth_mapq,
        exclude_dup => $exclude_dup,
    );

    my ($gene_mean, $gene_len, $win_depth_ref) = calculate_depth_from_gene_depth_file(
        depth_file  => $depth_file,
        gene        => $gene,
        windows_ref => $gene_windows{$gkey},
    );

    $gene_mean_depth{$gkey}   = $gene_mean;
    $gene_total_length{$gkey} = $gene_len;

    foreach my $wkey (keys %$win_depth_ref) {
        $window_mean_depth{$wkey} = $win_depth_ref->{$wkey};
    }

    unlink $depth_file if !$keep_depth_file;
}

if (!$keep_depth_file) {
    rmdir $depth_dir if -d $depth_dir;
}

# ============================================================
# Step 4. Write gene mean depth output
# ============================================================

open my $GD, ">", $gene_depth_out
    or die "[ERROR] Cannot write gene mean depth output: $gene_depth_out\n";

print $GD join("\t", qw(
    SampleID Gene Classification Transcript Chrom Start End Strand ExonCount
    Gene_Key Gene_Length Gene_Mean_Depth Window_Count
)), "\n";

foreach my $gene (@genes) {
    my $gkey = gene_key($gene);
    my $window_count = exists $gene_windows{$gkey}
        ? scalar(@{ $gene_windows{$gkey} })
        : 0;

    print $GD join("\t",
        $sample,
        $gene->{gene},
        $gene->{classification},
        $gene->{transcript},
        $gene->{chr},
        $gene->{start},
        $gene->{end},
        $gene->{strand},
        $gene->{exon_count},
        $gkey,
        $gene_total_length{$gkey},
        sprintf("%.4f", $gene_mean_depth{$gkey}),
        $window_count,
    ), "\n";
}

close $GD;

# ============================================================
# Step 5. Write window depth output
# ============================================================

open my $WD, ">", $window_depth_out
    or die "[ERROR] Cannot write window depth output: $window_depth_out\n";

print $WD join("\t", qw(
    SampleID Gene Classification Transcript Window_ID Chrom Start End Length
    Gene_Start Gene_End Window_Mean_Depth
)), "\n";

foreach my $win (@windows) {
    my $wkey = window_key($win);

    print $WD join("\t",
        $sample,
        $win->{gene},
        $win->{classification},
        $win->{transcript},
        $win->{window_id},
        $win->{chr},
        $win->{start},
        $win->{end},
        $win->{length},
        $win->{gene_start},
        $win->{gene_end},
        sprintf("%.4f", $window_mean_depth{$wkey} // 0),
    ), "\n";
}

close $WD;

# ============================================================
# Step 6. Calculate ratio and classify windows
# ============================================================

my @ratio_records;

open my $RATIO, ">", $all_ratio_out
    or die "[ERROR] Cannot write all window ratio output: $all_ratio_out\n";

print $RATIO join("\t", qw(
    SampleID Gene Classification Transcript Window_ID Chrom Start End Length
    Gene_Start Gene_End Gene_Mean_Depth Window_Mean_Depth
    Depth_Ratio Window_Status Comment
)), "\n";

foreach my $win (@windows) {
    my $wkey = window_key($win);
    my $gkey = $win->{gene_key};

    my $gene_mean   = $gene_mean_depth{$gkey};
    my $window_mean = $window_mean_depth{$wkey} // 0;

    my $ratio = "NA";
    $ratio = $window_mean / $gene_mean if defined $gene_mean && $gene_mean > 0;

    my ($status, $comment) = classify_window_depth(
        gene_mean_depth       => $gene_mean,
        window_mean_depth     => $window_mean,
        ratio                 => $ratio,
        min_gene_mean_depth   => $min_gene_mean_depth,
        min_window_mean_depth => $min_window_mean_depth,
        depth_ratio_cutoff    => $depth_ratio_cutoff,
        het_del_ratio_low     => $het_del_ratio_low,
        het_del_ratio_high    => $het_del_ratio_high,
    );

    my $record = {
        %$win,
        gene_mean_depth   => $gene_mean,
        window_mean_depth => $window_mean,
        depth_ratio       => $ratio,
        window_status     => $status,
        comment           => $comment,
    };

    push @ratio_records, $record;

    print $RATIO join("\t",
        $sample,
        $win->{gene},
        $win->{classification},
        $win->{transcript},
        $win->{window_id},
        $win->{chr},
        $win->{start},
        $win->{end},
        $win->{length},
        $win->{gene_start},
        $win->{gene_end},
        sprintf("%.4f", $gene_mean),
        sprintf("%.4f", $window_mean),
        ($ratio eq "NA" ? "NA" : sprintf("%.4f", $ratio)),
        $status,
        $comment,
    ), "\n";
}

close $RATIO;

# ============================================================
# Step 7. Merge low-depth windows into candidates
# ============================================================

my @candidates = merge_low_depth_windows(
    records               => \@ratio_records,
    min_candidate_windows => $min_candidate_windows,
    max_merge_gap_windows => $max_merge_gap_windows,
    min_candidate_size    => $min_depth_candidate_size,
    max_candidate_size    => $max_depth_candidate_size,
);

# ============================================================
# Step 8. Write candidate output
# ============================================================

open my $CAND, ">", $out
    or die "[ERROR] Cannot write candidate output: $out\n";

print $CAND join("\t", qw(
    SampleID Gene Classification Transcript Chrom Start End Candidate_Length
    Window_Count Gene_Mean_Depth Candidate_Mean_Depth
    Mean_Depth_Ratio Min_Depth_Ratio Max_Depth_Ratio
    Depth_Status Candidate_Status Comment
)), "\n";

foreach my $cand (@candidates) {
    print $CAND join("\t",
        $sample,
        $cand->{gene},
        $cand->{classification},
        $cand->{transcript},
        $cand->{chr},
        $cand->{start},
        $cand->{end},
        $cand->{candidate_length},
        $cand->{window_count},
        sprintf("%.4f", $cand->{gene_mean_depth}),
        sprintf("%.4f", $cand->{candidate_mean_depth}),
        sprintf("%.4f", $cand->{mean_depth_ratio}),
        sprintf("%.4f", $cand->{min_depth_ratio}),
        sprintf("%.4f", $cand->{max_depth_ratio}),
        $cand->{depth_status},
        "Depth_candidate",
        $cand->{comment},
    ), "\n";
}

close $CAND;

print "[INFO] Window-based depth analysis finished\n";
print "[INFO] Sample                  : $sample\n";
print "[INFO] Gene TXT                : $gene_txt\n";
print "[INFO] Window depth output     : $window_depth_out\n";
print "[INFO] Gene mean depth output  : $gene_depth_out\n";
print "[INFO] All window ratio output : $all_ratio_out\n";
print "[INFO] Candidate output        : $out\n";
print "[INFO] Candidate number        : " . scalar(@candidates) . "\n";

exit 0;

# ============================================================
# Subroutines
# ============================================================

sub read_config {
    my ($file) = @_;
    my %conf;

    open my $FH, "<", $file or die "[ERROR] Cannot open config file: $file\n";

    while (my $line = <$FH>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        if ($line =~ /^\s*([^=\s]+)\s*=\s*(.*?)\s*$/) {
            my $key = $1;
            my $val = $2;

            $val =~ s/\s+#.*$//;
            $val =~ s/^\s+|\s+$//g;
            $val =~ s/^['"]//;
            $val =~ s/['"]$//;

            $conf{$key} = $val;
        }
    }

    close $FH;
    return %conf;
}

sub get_conf {
    my ($conf_ref, $key, $default) = @_;

    if (exists $conf_ref->{$key} && defined $conf_ref->{$key} && $conf_ref->{$key} ne "") {
        return $conf_ref->{$key};
    }

    return $default;
}

sub resolve_path_by_config {
    my ($path, $config_file) = @_;

    die "[ERROR] Empty path provided\n"
        unless defined $path && $path ne "";

    if ($path =~ /^\//) {
        my $abs = abs_path($path);
        die "[ERROR] Path not found: $path\n" unless defined $abs;
        return $abs;
    }

    my $conf_dir = dirname($config_file);

    my $try1 = "$conf_dir/$path";
    my $abs1 = abs_path($try1);
    return $abs1 if defined $abs1;

    my $project_dir = abs_path("$conf_dir/..");
    my $try2 = "$project_dir/$path";
    my $abs2 = abs_path($try2);
    return $abs2 if defined $abs2;

    die "[ERROR] Path not found: $path\n"
      . "        Tried: $try1\n"
      . "        Tried: $try2\n";
}

sub read_core_gene_list {
    my ($file) = @_;

    my %gene;
    my %classification;

    open my $FH, "<", $file or die "[ERROR] Cannot open HCM core gene list: $file\n";

    my $line_no = 0;

    while (my $line = <$FH>) {
        chomp $line;
        $line =~ s/\r$//;
        $line =~ s/^\s+|\s+$//g;
        $line_no++;

        next if $line eq "";
        next if $line =~ /^#/;

        my @F = split /\t|\s+/, $line;
        my $g = $F[0] // "";
        my $c = $F[1] // "NA";

        next if $line_no == 1 && $g =~ /^Gene$/i;

        next unless $g ne "";

        $gene{$g} = 1;

        if (!exists $classification{$g}) {
            $classification{$g} = $c;
        } else {
            $classification{$g} = choose_higher_classification($classification{$g}, $c);
        }
    }

    close $FH;

    die "[ERROR] No genes loaded from HCM core gene list: $file\n"
        unless keys %gene;

    return (%gene, %classification);
}

sub choose_higher_classification {
    my ($old, $new) = @_;

    my %rank = (
        "Definitive" => 5,
        "Strong"     => 4,
        "Moderate"   => 3,
        "Limited"    => 2,
        "Disputed"   => 1,
        "NA"         => 0,
    );

    $old = "NA" unless defined $old && $old ne "";
    $new = "NA" unless defined $new && $new ne "";

    my $old_rank = exists $rank{$old} ? $rank{$old} : 0;
    my $new_rank = exists $rank{$new} ? $rank{$new} : 0;

    return $new_rank > $old_rank ? $new : $old;
}

sub read_gene_txt {
    my %args = @_;

    my $file               = $args{file};
    my $analyze_core_only  = $args{analyze_core_only} || 0;
    my $core_gene_ref      = $args{core_gene_ref} || {};
    my $classification_ref = $args{classification_ref} || {};

    my @genes;
    my %col;

    open my $FH, "<", $file or die "[ERROR] Cannot open gene TXT: $file\n";

    my $line_no = 0;

    while (my $line = <$FH>) {
        chomp $line;
        $line =~ s/\r$//;
        $line_no++;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        my @F = split /\t/, $line;

        if ($line_no == 1) {
            for my $i (0 .. $#F) {
                $col{$F[$i]} = $i;
            }

            foreach my $required (qw(Gene Transcript Chrom Start End Strand ExonCount)) {
                die "[ERROR] Gene TXT missing required column '$required': $file\n"
                    unless exists $col{$required};
            }

            next;
        }

        my $gene       = $F[$col{Gene}];
        my $transcript = $F[$col{Transcript}];
        my $chr        = $F[$col{Chrom}];
        my $start      = $F[$col{Start}];
        my $end        = $F[$col{End}];
        my $strand     = $F[$col{Strand}];
        my $exon_count = $F[$col{ExonCount}];

        next unless defined $gene && defined $chr && defined $start && defined $end;

        if ($analyze_core_only) {
            next unless exists $core_gene_ref->{$gene};
        }

        next unless $start =~ /^\d+$/ && $end =~ /^\d+$/;
        next if $end < $start;

        my $classification = exists $classification_ref->{$gene}
            ? $classification_ref->{$gene}
            : "NA";

        push @genes, {
            gene           => $gene,
            classification => $classification,
            transcript     => $transcript,
            chr            => $chr,
            start          => $start,
            end            => $end,
            strand         => $strand,
            exon_count     => $exon_count,
            length         => $end - $start + 1,
        };
    }

    close $FH;
    return @genes;
}

sub generate_windows {
    my %args = @_;

    my $gene            = $args{gene};
    my $sample          = $args{sample};
    my $window_size     = $args{window_size};
    my $window_step     = $args{window_step};
    my $min_window_size = $args{min_window_size};

    my @windows;
    my $gkey = gene_key($gene);
    my $window_index = 0;

    for (my $start = $gene->{start}; $start <= $gene->{end}; $start += $window_step) {
        my $end = $start + $window_size - 1;
        $end = $gene->{end} if $end > $gene->{end};

        my $length = $end - $start + 1;
        next if $length < $min_window_size;

        $window_index++;

        push @windows, {
            sample         => $sample,
            gene           => $gene->{gene},
            classification => $gene->{classification},
            transcript     => $gene->{transcript},
            chr            => $gene->{chr},

            start          => $start,
            end            => $end,
            length         => $length,

            gene_start     => $gene->{start},
            gene_end       => $gene->{end},
            gene_length    => $gene->{length},

            strand         => $gene->{strand},
            exon_count     => $gene->{exon_count},

            gene_key       => $gkey,
            window_id      => $gene->{gene} . "_W" . $window_index,
        };
    }

    return @windows;
}

sub run_samtools_depth_for_gene {
    my %args = @_;

    my $samtools    = $args{samtools};
    my $bam         = $args{bam};
    my $gene        = $args{gene};
    my $depth_file  = $args{depth_file};
    my $min_mapq    = $args{min_mapq};
    my $exclude_dup = $args{exclude_dup};

    my $region = join_region($gene->{chr}, $gene->{start}, $gene->{end});

    my @cmd = (
        $samtools,
        "depth",
        "-a",
        "-q", $min_mapq,
        "-r", $region,
    );

    if ($exclude_dup) {
        push @cmd, ("-G", "1024");
    }

    push @cmd, $bam;

    my $cmd = join(" ", @cmd) . " > $depth_file";

    system($cmd) == 0
        or die "[ERROR] Failed to run command:\n$cmd\n";
}

sub calculate_depth_from_gene_depth_file {
    my %args = @_;

    my $depth_file  = $args{depth_file};
    my $gene        = $args{gene};
    my $windows_ref = $args{windows_ref};

    my @windows = @$windows_ref;

    my %window_depth_sum;
    my %window_length;

    foreach my $win (@windows) {
        my $wkey = window_key($win);
        $window_depth_sum{$wkey} = 0;
        $window_length{$wkey}    = $win->{length};
    }

    my $gene_depth_sum = 0;
    my $gene_length = $gene->{end} - $gene->{start} + 1;

    open my $IN, "<", $depth_file
        or die "[ERROR] Cannot open depth file: $depth_file\n";

    while (my $line = <$IN>) {
        chomp $line;
        next if $line =~ /^\s*$/;

        my @F = split /\t/, $line;
        next unless @F >= 3;

        my ($chr, $pos, $depth) = @F[0, 1, 2];

        next unless defined $pos && $pos =~ /^\d+$/;
        next unless defined $depth && $depth =~ /^\d+$/;

        $gene_depth_sum += $depth;

        foreach my $win (@windows) {
            next if $pos < $win->{start};
            last if $pos > $win->{end} && $win->{start} > $pos;
            next if $pos > $win->{end};

            my $wkey = window_key($win);
            $window_depth_sum{$wkey} += $depth;
        }
    }

    close $IN;

    my $gene_mean_depth = 0;
    $gene_mean_depth = $gene_depth_sum / $gene_length if $gene_length > 0;

    my %window_mean_depth;

    foreach my $win (@windows) {
        my $wkey = window_key($win);
        my $mean = 0;

        if ($window_length{$wkey} && $window_length{$wkey} > 0) {
            $mean = $window_depth_sum{$wkey} / $window_length{$wkey};
        }

        $window_mean_depth{$wkey} = $mean;
    }

    return ($gene_mean_depth, $gene_length, \%window_mean_depth);
}

sub classify_window_depth {
    my %args = @_;

    my $gene_mean_depth       = $args{gene_mean_depth};
    my $window_mean_depth     = $args{window_mean_depth};
    my $ratio                 = $args{ratio};
    my $min_gene_mean_depth   = $args{min_gene_mean_depth};
    my $min_window_mean_depth = $args{min_window_mean_depth};
    my $cutoff                = $args{depth_ratio_cutoff};
    my $het_low               = $args{het_del_ratio_low};
    my $het_high              = $args{het_del_ratio_high};

    if (!defined $gene_mean_depth || $gene_mean_depth < $min_gene_mean_depth) {
        return ("Low_gene_depth", "Gene mean depth is lower than threshold");
    }

    if (!defined $window_mean_depth || $window_mean_depth < $min_window_mean_depth) {
        return ("Low_absolute_window_depth", "Window mean depth is lower than minimum evaluable depth");
    }

    if (!defined $ratio || $ratio eq "NA") {
        return ("NA", "Depth ratio cannot be calculated");
    }

    if ($ratio <= $cutoff) {
        if ($ratio >= $het_low && $ratio <= $het_high) {
            return ("Low_depth_window", "Depth ratio is consistent with possible heterozygous deletion");
        } elsif ($ratio < $het_low) {
            return ("Very_low_depth_window", "Depth ratio is lower than expected heterozygous deletion range");
        } else {
            return ("Mild_low_depth_window", "Depth ratio is below cutoff");
        }
    }

    return ("Normal_depth_window", "Depth ratio is not reduced");
}

sub merge_low_depth_windows {
    my %args = @_;

    my @records = @{ $args{records} };
    my $min_windows = $args{min_candidate_windows};
    my $max_gap = $args{max_merge_gap_windows};
    my $min_size = $args{min_candidate_size};
    my $max_size = $args{max_candidate_size};

    my %grouped;

    foreach my $r (@records) {
        my $key = join("|",
            $r->{gene},
            $r->{transcript},
            $r->{chr},
            $r->{gene_start},
            $r->{gene_end},
        );

        push @{ $grouped{$key} }, $r;
    }

    my @candidates;

    foreach my $gkey (sort keys %grouped) {
        my @sorted = sort { $a->{start} <=> $b->{start} } @{ $grouped{$gkey} };

        my @current;
        my $gap_count = 0;

        foreach my $r (@sorted) {
            my $is_low = is_low_depth_status($r->{window_status});

            if ($is_low) {
                push @current, $r;
                $gap_count = 0;
            } else {
                if (@current && $gap_count < $max_gap) {
                    push @current, $r;
                    $gap_count++;
                } else {
                    push_candidate_if_valid(
                        candidates_ref => \@candidates,
                        windows_ref    => \@current,
                        min_windows    => $min_windows,
                        min_size       => $min_size,
                        max_size       => $max_size,
                    ) if @current;

                    @current = ();
                    $gap_count = 0;
                }
            }
        }

        push_candidate_if_valid(
            candidates_ref => \@candidates,
            windows_ref    => \@current,
            min_windows    => $min_windows,
            min_size       => $min_size,
            max_size       => $max_size,
        ) if @current;
    }

    return @candidates;
}

sub push_candidate_if_valid {
    my %args = @_;

    my $candidates_ref = $args{candidates_ref};
    my $windows_ref    = $args{windows_ref};
    my $min_windows    = $args{min_windows};
    my $min_size       = $args{min_size};
    my $max_size       = $args{max_size};

    my @low_windows = grep { is_low_depth_status($_->{window_status}) } @$windows_ref;
    return if scalar(@low_windows) < $min_windows;

    my $candidate = build_candidate(\@low_windows);

    return if $candidate->{candidate_length} < $min_size;
    return if $candidate->{candidate_length} > $max_size;

    push @$candidates_ref, $candidate;
}

sub is_low_depth_status {
    my ($status) = @_;

    return 1 if $status eq "Low_depth_window";
    return 1 if $status eq "Very_low_depth_window";
    return 1 if $status eq "Mild_low_depth_window";

    return 0;
}

sub build_candidate {
    my ($wins_ref) = @_;

    my @wins = sort { $a->{start} <=> $b->{start} } @$wins_ref;

    my $first = $wins[0];
    my $last  = $wins[-1];

    my $start = $first->{start};
    my $end   = $last->{end};

    my $total_depth_sum = 0;
    my $total_length = 0;

    my $ratio_sum = 0;
    my $ratio_count = 0;
    my $min_ratio = 999999;
    my $max_ratio = -1;

    foreach my $w (@wins) {
        $total_depth_sum += $w->{window_mean_depth} * $w->{length};
        $total_length    += $w->{length};

        if (defined $w->{depth_ratio} && $w->{depth_ratio} ne "NA") {
            $ratio_sum += $w->{depth_ratio};
            $ratio_count++;

            $min_ratio = $w->{depth_ratio} if $w->{depth_ratio} < $min_ratio;
            $max_ratio = $w->{depth_ratio} if $w->{depth_ratio} > $max_ratio;
        }
    }

    my $candidate_mean_depth = 0;
    $candidate_mean_depth = $total_depth_sum / $total_length if $total_length > 0;

    my $mean_depth_ratio = 0;
    $mean_depth_ratio = $ratio_sum / $ratio_count if $ratio_count > 0;

    my $depth_status = "Low_depth_region";
    my $comment = "Merged low-depth windows";

    if ($mean_depth_ratio >= 0.35 && $mean_depth_ratio <= 0.70) {
        $depth_status = "Heterozygous_deletion_like";
        $comment = "Merged low-depth windows are consistent with possible heterozygous deletion";
    } elsif ($mean_depth_ratio < 0.35) {
        $depth_status = "Severe_depth_reduction";
        $comment = "Merged low-depth windows show severe depth reduction";
    }

    return {
        gene                 => $first->{gene},
        classification       => $first->{classification},
        transcript           => $first->{transcript},
        chr                  => $first->{chr},
        start                => $start,
        end                  => $end,
        candidate_length     => $end - $start + 1,
        window_count         => scalar(@wins),
        gene_mean_depth      => $first->{gene_mean_depth},
        candidate_mean_depth => $candidate_mean_depth,
        mean_depth_ratio     => $mean_depth_ratio,
        min_depth_ratio      => $min_ratio == 999999 ? 0 : $min_ratio,
        max_depth_ratio      => $max_ratio == -1 ? 0 : $max_ratio,
        depth_status         => $depth_status,
        comment              => $comment,
    };
}

sub make_gene_depth_filename {
    my %args = @_;

    my $depth_dir = $args{depth_dir};
    my $sample    = $args{sample};
    my $gene      = $args{gene};

    my $name = join("_",
        $sample,
        $gene->{gene},
        $gene->{transcript},
        $gene->{chr},
        $gene->{start},
        $gene->{end},
    );

    $name =~ s/[^\w.\-]+/_/g;

    return "$depth_dir/$name.depth.tsv";
}

sub gene_key {
    my ($r) = @_;

    my $start = exists $r->{gene_start} ? $r->{gene_start} : $r->{start};
    my $end   = exists $r->{gene_end}   ? $r->{gene_end}   : $r->{end};

    return join("|",
        $r->{gene},
        $r->{transcript},
        $r->{chr},
        $start,
        $end,
    );
}

sub window_key {
    my ($win) = @_;

    return join("|",
        $win->{gene},
        $win->{transcript},
        $win->{window_id},
        $win->{chr},
        $win->{start},
        $win->{end},
    );
}

sub join_region {
    my ($chr, $start, $end) = @_;
    return $chr . ":" . $start . "-" . $end;
}

sub check_positive_integer {
    my ($name, $value) = @_;

    die "[ERROR] $name must be a positive integer\n"
        unless defined $value && $value =~ /^\d+$/ && $value > 0;
}

sub check_non_negative_integer {
    my ($name, $value) = @_;

    die "[ERROR] $name must be a non-negative integer\n"
        unless defined $value && $value =~ /^\d+$/;
}

sub usage {
    return <<"USAGE";

Usage:
  perl bin/run_gene_mean_depth.pl \\
      --config conf/hcm_exondel.example.conf \\
      --bam sample.sorted.bam \\
      --sample SAMPLE001 \\
      --out results/SAMPLE001/01.depth/SAMPLE001.depth_candidates.tsv \\
      --threads 4

Required arguments:
  --config     Config file.
  --bam        Coordinate-sorted BAM file.
  --sample     Sample ID.
  --out        Output depth candidate file.

Optional arguments:
  --threads    Number of threads. Reserved for compatibility.
  --help       Show this help message.

Required config:
  SAMTOOLS
  REFSEQ_MANE_SELECT_GENE_TXT

Recommended config:
  HCM_CORE_GENE_LIST
  ANALYZE_CORE_GENES_ONLY=1

Core gene list format:
  Gene    Classification
  FHOD3   Definitive
  MYH7    Definitive

Gene TXT format:
  Gene    Transcript    Chrom    Start    End    Strand    ExonCount

Coordinate:
  1-based closed interval.

Main output:
  depth_candidates.tsv

Additional output:
  *.gene_mean_depth.tsv
  *.window_depth.tsv
  *.all_window_ratio.tsv

Intermediate output:
  *.gene_depth_files/*.depth.tsv

USAGE
}


