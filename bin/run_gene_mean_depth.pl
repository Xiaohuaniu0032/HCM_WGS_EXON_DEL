#!/usr/bin/env perl
# -*- coding: utf-8 -*-

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path);
use IO::Handle;

STDOUT->autoflush(1);
STDERR->autoflush(1);

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
#   5. If a per-gene depth file already exists and is non-empty,
#      it will be reused automatically.
#   6. Gene mean depth and window mean depth are calculated
#      from per-gene base-level depth file.
#   7. Window depth is calculated using prefix sum for speed.
#   8. A window is considered deletion-supporting when:
#        Depth_Ratio <= DEL_DEPTH_RATIO_CUTOFF
#   9. A final depth candidate is reported only when at least:
#        MIN_CONSECUTIVE_DEL_WINDOWS
#      consecutive deletion-supporting windows are observed.
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
#
# Recommended config:
#   ANALYZE_CORE_GENES_ONLY=1
#   HCM_CORE_GENE_LIST
#   WINDOW_SIZE=1000
#   WINDOW_STEP=500
#   MIN_WINDOW_SIZE=100
#   MIN_GENE_MEAN_DEPTH=20
#   MIN_DEPTH_MAPQ=20
#   EXCLUDE_DUPLICATES_FOR_DEPTH=1
#   DEL_DEPTH_RATIO_CUTOFF=0.65
#   MIN_CONSECUTIVE_DEL_WINDOWS=3
#   KEEP_GENE_DEPTH_FILE=1
#
# Output:
#   1. *.gene_mean_depth.tsv
#   2. *.window_depth.tsv
#   3. *.all_window_ratio.tsv
#   4. *.del_windows.tsv
#   5. *.depth_candidates.tsv
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

    my ($core_gene_ref, $gene_classification_ref) = read_core_gene_list($core_gene_list);

    %core_gene           = %{$core_gene_ref};
    %gene_classification = %{$gene_classification_ref};
}

my $window_size     = get_conf(\%CONF, "WINDOW_SIZE", 1000);
my $window_step     = get_conf(\%CONF, "WINDOW_STEP", 500);
my $min_window_size = get_conf(\%CONF, "MIN_WINDOW_SIZE", 100);

my $min_gene_mean_depth = get_conf(\%CONF, "MIN_GENE_MEAN_DEPTH", 20);

my $min_depth_mapq = get_conf(\%CONF, "MIN_DEPTH_MAPQ", get_conf(\%CONF, "MIN_MAPQ", 20));
my $exclude_dup    = get_conf(\%CONF, "EXCLUDE_DUPLICATES_FOR_DEPTH", 1);

my $del_depth_ratio_cutoff = get_conf(\%CONF, "DEL_DEPTH_RATIO_CUTOFF", 0.65);
my $min_consecutive_del_windows = get_conf(\%CONF, "MIN_CONSECUTIVE_DEL_WINDOWS", 3);

# Default is 1 because existing gene depth files are reused automatically.
my $keep_depth_file = get_conf(\%CONF, "KEEP_GENE_DEPTH_FILE", 1);

check_positive_integer("WINDOW_SIZE", $window_size);
check_positive_integer("WINDOW_STEP", $window_step);
check_positive_integer("MIN_WINDOW_SIZE", $min_window_size);
check_non_negative_number("MIN_GENE_MEAN_DEPTH", $min_gene_mean_depth);
check_non_negative_integer("MIN_DEPTH_MAPQ", $min_depth_mapq);
check_binary_flag("EXCLUDE_DUPLICATES_FOR_DEPTH", $exclude_dup);
check_numeric_range("DEL_DEPTH_RATIO_CUTOFF", $del_depth_ratio_cutoff, 0, 1);
check_positive_integer("MIN_CONSECUTIVE_DEL_WINDOWS", $min_consecutive_del_windows);
check_binary_flag("KEEP_GENE_DEPTH_FILE", $keep_depth_file);

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
my $del_window_out   = $prefix . ".del_windows.tsv";

# ============================================================
# Step 1. Read gene regions
# ============================================================

my @genes = read_gene_txt(
    file               => $gene_txt,
    analyze_core_only  => $analyze_core_only,
    core_gene_ref      => \%core_gene,
    classification_ref => \%gene_classification,
);

die "[ERROR] No valid gene regions found in: $gene_txt\n" unless @genes;

print "[INFO] Window-based depth analysis started\n";
print "[INFO] Sample                         : $sample\n";
print "[INFO] BAM                            : $bam\n";
print "[INFO] Config                         : $config\n";
print "[INFO] Gene coordinate file           : $gene_txt\n";
print "[INFO] WINDOW_SIZE                    : $window_size\n";
print "[INFO] WINDOW_STEP                    : $window_step\n";
print "[INFO] MIN_WINDOW_SIZE                : $min_window_size\n";
print "[INFO] MIN_GENE_MEAN_DEPTH            : $min_gene_mean_depth\n";
print "[INFO] MIN_DEPTH_MAPQ                 : $min_depth_mapq\n";
print "[INFO] EXCLUDE_DUPLICATES_FOR_DEPTH   : $exclude_dup\n";
print "[INFO] DEL_DEPTH_RATIO_CUTOFF         : $del_depth_ratio_cutoff\n";
print "[INFO] MIN_CONSECUTIVE_DEL_WINDOWS    : $min_consecutive_del_windows\n";
print "[INFO] KEEP_GENE_DEPTH_FILE           : $keep_depth_file\n";

if ($analyze_core_only) {
    print "[INFO] ANALYZE_CORE_GENES_ONLY        : 1\n";
    print "[INFO] Core gene list                 : $core_gene_list\n";
    print "[INFO] Core genes loaded              : " . scalar(keys %core_gene) . "\n";
    print "[INFO] Gene classifications loaded    : " . scalar(keys %gene_classification) . "\n";
    print "[INFO] Gene regions retained          : " . scalar(@genes) . "\n";
} else {
    print "[INFO] ANALYZE_CORE_GENES_ONLY        : 0\n";
    print "[INFO] Gene regions retained          : " . scalar(@genes) . "\n";
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

print "[INFO] Total windows generated        : " . scalar(@windows) . "\n";

# ============================================================
# Step 3. Generate or reuse per-gene depth file and calculate depths
# ============================================================

my %gene_mean_depth;
my %gene_total_length;
my %window_mean_depth;

my $reused_depth_file_count = 0;
my $generated_depth_file_count = 0;

my $gene_index = 0;
my $gene_total = scalar(@genes);

foreach my $gene (@genes) {
    $gene_index++;

    my $gkey = gene_key($gene);

    my $depth_file = make_gene_depth_filename(
        depth_dir => $depth_dir,
        sample    => $sample,
        gene      => $gene,
    );

    print "[INFO] [$gene_index/$gene_total] Processing gene: "
        . $gene->{gene} . " "
        . $gene->{transcript} . " "
        . $gene->{chr} . ":"
        . $gene->{start} . "-"
        . $gene->{end} . "\n";

    if (-s $depth_file) {
        print "[INFO] [$gene_index/$gene_total] Reuse existing gene depth file: $depth_file\n";
        $reused_depth_file_count++;
    } else {
        print "[INFO] [$gene_index/$gene_total] Gene depth file not found or empty, generate now: $depth_file\n";

        run_samtools_depth_for_gene(
            samtools    => $samtools,
            bam         => $bam,
            gene        => $gene,
            depth_file  => $depth_file,
            min_mapq    => $min_depth_mapq,
            exclude_dup => $exclude_dup,
        );

        $generated_depth_file_count++;
    }

    print "[INFO] [$gene_index/$gene_total] Calculating gene/window depth from: $depth_file\n";

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

    print "[INFO] [$gene_index/$gene_total] Finished gene: "
        . $gene->{gene}
        . ", Gene_Mean_Depth="
        . sprintf("%.4f", $gene_mean)
        . "\n";

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
my $del_window_count = 0;

open my $RATIO, ">", $all_ratio_out
    or die "[ERROR] Cannot write all window ratio output: $all_ratio_out\n";

open my $DELWIN, ">", $del_window_out
    or die "[ERROR] Cannot write deletion window output: $del_window_out\n";

my @ratio_header = qw(
    SampleID Gene Classification Transcript Window_ID Chrom Start End Length
    Gene_Start Gene_End Gene_Mean_Depth Window_Mean_Depth
    Depth_Ratio Is_Del_Window Window_Status Comment
);

print $RATIO  join("\t", @ratio_header), "\n";
print $DELWIN join("\t", @ratio_header), "\n";

foreach my $win (@windows) {
    my $wkey = window_key($win);
    my $gkey = $win->{gene_key};

    my $gene_mean   = $gene_mean_depth{$gkey};
    my $window_mean = $window_mean_depth{$wkey} // 0;

    my $ratio = "NA";
    $ratio = $window_mean / $gene_mean
        if defined $gene_mean && $gene_mean > 0;

    my ($status, $comment, $is_del_window) = classify_window_depth(
        gene_mean_depth     => $gene_mean,
        ratio               => $ratio,
        min_gene_mean_depth => $min_gene_mean_depth,
        del_ratio_cutoff    => $del_depth_ratio_cutoff,
    );

    my $record = {
        %$win,
        gene_mean_depth   => $gene_mean,
        window_mean_depth => $window_mean,
        depth_ratio       => $ratio,
        is_del_window     => $is_del_window,
        window_status     => $status,
        comment           => $comment,
    };

    push @ratio_records, $record;

    my @out_fields = (
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
        $is_del_window ? "Yes" : "No",
        $status,
        $comment,
    );

    print $RATIO join("\t", @out_fields), "\n";

    if ($is_del_window) {
        print $DELWIN join("\t", @out_fields), "\n";
        $del_window_count++;
    }
}

close $RATIO;
close $DELWIN;

# ============================================================
# Step 7. Merge consecutive del windows into candidates
# ============================================================

my @candidates = merge_consecutive_del_windows(
    records                     => \@ratio_records,
    min_consecutive_del_windows => $min_consecutive_del_windows,
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
    Candidate_Status Del_Window_IDs Comment
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
        "Depth_candidate",
        $cand->{del_window_ids},
        $cand->{comment},
    ), "\n";
}

close $CAND;

print "[INFO] Window-based depth analysis finished\n";
print "[INFO] Sample                         : $sample\n";
print "[INFO] Gene mean depth output         : $gene_depth_out\n";
print "[INFO] Window depth output            : $window_depth_out\n";
print "[INFO] All window ratio output        : $all_ratio_out\n";
print "[INFO] Del window output              : $del_window_out\n";
print "[INFO] Candidate output               : $out\n";
print "[INFO] Gene depth files reused        : $reused_depth_file_count\n";
print "[INFO] Gene depth files generated     : $generated_depth_file_count\n";
print "[INFO] Del window number              : $del_window_count\n";
print "[INFO] Candidate number               : " . scalar(@candidates) . "\n";

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

        $g =~ s/^\s+|\s+$//g;
        $c =~ s/^\s+|\s+$//g;
        $c = "NA" if $c eq "";

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

    return (\%gene, \%classification);
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

        $gene       =~ s/^\s+|\s+$//g;
        $transcript =~ s/^\s+|\s+$//g if defined $transcript;
        $chr        =~ s/^\s+|\s+$//g;
        $strand     =~ s/^\s+|\s+$//g if defined $strand;

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
        shell_quote($samtools),
        "depth",
        "-a",
        "-q", shell_quote($min_mapq),
        "-r", shell_quote($region),
    );

    if ($exclude_dup) {
        push @cmd, ("-G", "1024");
    }

    push @cmd, shell_quote($bam);

    my $cmd = join(" ", @cmd) . " > " . shell_quote($depth_file);

    print "[INFO] Generate gene depth file: $depth_file\n";
    print "[INFO] Command: $cmd\n";

    system($cmd) == 0
        or die "[ERROR] Failed to run command:\n$cmd\n";
}

sub calculate_depth_from_gene_depth_file {
    my %args = @_;

    my $depth_file  = $args{depth_file};
    my $gene        = $args{gene};
    my $windows_ref = $args{windows_ref};

    my @windows = @$windows_ref;

    my $gene_start  = $gene->{start};
    my $gene_end    = $gene->{end};
    my $gene_length = $gene_end - $gene_start + 1;

    die "[ERROR] Invalid gene length for "
        . $gene->{gene} . ": $gene_start-$gene_end\n"
        if $gene_length <= 0;

    # depth_by_offset[1] corresponds to gene_start.
    # depth_by_offset[2] corresponds to gene_start + 1.
    my @depth_by_offset;

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

        next if $pos < $gene_start;
        next if $pos > $gene_end;

        my $offset = $pos - $gene_start + 1;
        $depth_by_offset[$offset] = $depth;
    }

    close $IN;

    my @prefix_sum;
    $prefix_sum[0] = 0;

    for (my $i = 1; $i <= $gene_length; $i++) {
        my $d = defined $depth_by_offset[$i] ? $depth_by_offset[$i] : 0;
        $prefix_sum[$i] = $prefix_sum[$i - 1] + $d;
    }

    my $gene_mean_depth = 0;
    $gene_mean_depth = $prefix_sum[$gene_length] / $gene_length
        if $gene_length > 0;

    my %window_mean_depth;

    foreach my $win (@windows) {
        my $wkey = window_key($win);

        my $w_start = $win->{start};
        my $w_end   = $win->{end};

        $w_start = $gene_start if $w_start < $gene_start;
        $w_end   = $gene_end   if $w_end > $gene_end;

        my $w_len = $w_end - $w_start + 1;

        if ($w_len <= 0) {
            $window_mean_depth{$wkey} = 0;
            next;
        }

        my $left_offset  = $w_start - $gene_start + 1;
        my $right_offset = $w_end   - $gene_start + 1;

        my $window_depth_sum =
            $prefix_sum[$right_offset] - $prefix_sum[$left_offset - 1];

        my $mean = $window_depth_sum / $w_len;

        $window_mean_depth{$wkey} = $mean;
    }

    return ($gene_mean_depth, $gene_length, \%window_mean_depth);
}

sub classify_window_depth {
    my %args = @_;

    my $gene_mean_depth     = $args{gene_mean_depth};
    my $ratio               = $args{ratio};
    my $min_gene_mean_depth = $args{min_gene_mean_depth};
    my $cutoff              = $args{del_ratio_cutoff};

    if (!defined $gene_mean_depth || $gene_mean_depth < $min_gene_mean_depth) {
        return (
            "Low_gene_depth",
            "Gene mean depth is lower than MIN_GENE_MEAN_DEPTH; window is not evaluated as deletion",
            0
        );
    }

    if (!defined $ratio || $ratio eq "NA") {
        return (
            "NA",
            "Depth ratio cannot be calculated",
            0
        );
    }

    if ($ratio <= $cutoff) {
        return (
            "Del_window",
            "Depth ratio <= DEL_DEPTH_RATIO_CUTOFF",
            1
        );
    }

    return (
        "Normal_window",
        "Depth ratio > DEL_DEPTH_RATIO_CUTOFF",
        0
    );
}

sub merge_consecutive_del_windows {
    my %args = @_;

    my @records = @{ $args{records} };
    my $min_consecutive = $args{min_consecutive_del_windows};

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
        my @sorted = sort {
            $a->{start} <=> $b->{start}
                ||
            $a->{end} <=> $b->{end}
        } @{ $grouped{$gkey} };

        my @current_del_run;

        foreach my $r (@sorted) {
            if ($r->{is_del_window}) {
                push @current_del_run, $r;
            } else {
                add_del_run_candidate_if_valid(
                    candidates_ref  => \@candidates,
                    windows_ref     => \@current_del_run,
                    min_consecutive => $min_consecutive,
                );

                @current_del_run = ();
            }
        }

        add_del_run_candidate_if_valid(
            candidates_ref  => \@candidates,
            windows_ref     => \@current_del_run,
            min_consecutive => $min_consecutive,
        );
    }

    return @candidates;
}

sub add_del_run_candidate_if_valid {
    my %args = @_;

    my $candidates_ref  = $args{candidates_ref};
    my $windows_ref     = $args{windows_ref};
    my $min_consecutive = $args{min_consecutive};

    return unless @$windows_ref;
    return if scalar(@$windows_ref) < $min_consecutive;

    my $candidate = build_candidate($windows_ref);
    push @$candidates_ref, $candidate;
}

sub build_candidate {
    my ($wins_ref) = @_;

    my @wins = sort {
        $a->{start} <=> $b->{start}
            ||
        $a->{end} <=> $b->{end}
    } @$wins_ref;

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

    my @window_ids;

    foreach my $w (@wins) {
        push @window_ids, $w->{window_id};

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

    my $comment = "Merged consecutive deletion-supporting windows";

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
        del_window_ids       => join(",", @window_ids),
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

sub shell_quote {
    my ($s) = @_;
    die "[ERROR] Undefined value passed to shell_quote\n" unless defined $s;

    $s =~ s/'/'\\''/g;
    return "'$s'";
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

sub check_non_negative_number {
    my ($name, $value) = @_;

    die "[ERROR] $name must be a non-negative number\n"
        unless defined $value && $value =~ /^\d+(?:\.\d+)?$/ && $value >= 0;
}

sub check_numeric_range {
    my ($name, $value, $min, $max) = @_;

    die "[ERROR] $name must be numeric\n"
        unless defined $value && $value =~ /^\d+(?:\.\d+)?$/;

    die "[ERROR] $name must be >= $min and <= $max\n"
        if $value < $min || $value > $max;
}

sub check_binary_flag {
    my ($name, $value) = @_;

    die "[ERROR] $name must be 0 or 1\n"
        unless defined $value && $value =~ /^[01]$/;
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
  ANALYZE_CORE_GENES_ONLY=1
  HCM_CORE_GENE_LIST=db/hcm_core_genes.txt

Window config:
  WINDOW_SIZE=1000
  WINDOW_STEP=500
  MIN_WINDOW_SIZE=100

Depth QC config:
  MIN_GENE_MEAN_DEPTH=20
  MIN_DEPTH_MAPQ=20
  EXCLUDE_DUPLICATES_FOR_DEPTH=1

Deletion calling config:
  DEL_DEPTH_RATIO_CUTOFF=0.65
  MIN_CONSECUTIVE_DEL_WINDOWS=3

Intermediate file config:
  KEEP_GENE_DEPTH_FILE=1

Main output:
  *.depth_candidates.tsv

Additional output:
  *.gene_mean_depth.tsv
  *.window_depth.tsv
  *.all_window_ratio.tsv
  *.del_windows.tsv

Intermediate output:
  *.gene_depth_files/*.depth.tsv

Note:
  If *.gene_depth_files/*.depth.tsv already exists and is non-empty,
  this script will reuse it automatically and skip samtools depth.

  If MIN_DEPTH_MAPQ, EXCLUDE_DUPLICATES_FOR_DEPTH, BAM, or gene coordinates
  are changed, please delete the old *.gene_depth_files directory first.

USAGE
}


