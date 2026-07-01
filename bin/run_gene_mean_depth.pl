#!/usr/bin/env perl
# -*- coding: utf-8 -*-
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path);
use IO::Handle;

STDOUT->autoflush(1);
STDERR->autoflush(1);

# ============================================================
# run_gene_mean_depth.pl
#
# Purpose:
#   Window-based depth-ratio analysis for exon-level deletion
#   detection in HCM core genes.
#
# Pipeline logic:
#   1. HCMExonDel only analyzes HCM core genes.
#   2. This script always reads HCM_CORE_GENE_LIST and filters
#      REFSEQ_MANE_SELECT_GENE_TXT to those core genes.
#   3. The input BAM is expected to be SAMPLE.target.bam generated
#      by extract_target_bam.pl.
#   4. For each core gene, samtools depth is called once.
#   5. Gene mean depth and window mean depth are calculated from
#      per-gene base-level depth files.
#   6. Window depth is calculated using prefix sums.
#   7. A window is deletion-supporting when:
#        Depth_Ratio <= DEL_DEPTH_RATIO_CUTOFF
#   8. A depth candidate is reported when at least
#        MIN_CONSECUTIVE_DEL_WINDOWS
#      consecutive deletion-supporting windows are observed.
#
# Required command-line arguments:
#   --config
#   --bam
#   --sample
#   --out
#
# Required config parameters:
#   SAMTOOLS
#   REFSEQ_MANE_SELECT_GENE_TXT
#   HCM_CORE_GENE_LIST
#   WINDOW_SIZE
#   WINDOW_STEP
#   MIN_WINDOW_SIZE
#   MIN_GENE_MEAN_DEPTH
#   MIN_MAPQ
#   EXCLUDE_DUPLICATES
#   DEL_DEPTH_RATIO_CUTOFF
#   MIN_CONSECUTIVE_DEL_WINDOWS
#   KEEP_GENE_DEPTH_FILE
#
# Outputs:
#   1. *.gene_mean_depth.tsv
#   2. *.window_depth.tsv
#   3. *.all_window_ratio.tsv
#   4. *.del_windows.tsv
#   5. *.depth_candidates.tsv
#   6. *.gene_depth_files/*.depth.tsv, if KEEP_GENE_DEPTH_FILE=1
# ============================================================

my $config;
my $bam;
my $sample;
my $out;
my $help = 0;

GetOptions(
    "config=s" => \$config,
    "bam=s"    => \$bam,
    "sample=s" => \$sample,
    "out=s"    => \$out,
    "help"     => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage()
    unless defined $config
        && defined $bam
        && defined $sample
        && defined $out;

# ------------------------------------------------------------
# Validate runtime I/O
# ------------------------------------------------------------
$config = abs_path($config) if -e $config;
die "[ERROR] Config file not found: $config\n"
    unless defined $config && -s $config;
die "[ERROR] Config path must be absolute after resolution: $config\n"
    unless $config =~ m{^/};

$bam = abs_path($bam) if -e $bam;
die "[ERROR] BAM file not found: $bam\n"
    unless defined $bam && -s $bam;
die "[ERROR] BAM path must be absolute after resolution: $bam\n"
    unless $bam =~ m{^/};
die "[ERROR] BAM file must end with .bam: $bam\n"
    unless $bam =~ /\.bam$/i;

check_bam_index($bam);

die "[ERROR] Invalid sample name: $sample\n"
  . "        Sample name can only contain letters, numbers, dot, underscore and hyphen.\n"
    unless defined $sample && $sample =~ /^[A-Za-z0-9_.-]+$/;

die "[ERROR] Output path must be absolute: $out\n"
    unless defined $out && $out =~ m{^/};

my $outdir = dirname($out);
make_path($outdir) unless -d $outdir;
die "[ERROR] Cannot create output directory: $outdir\n"
    unless -d $outdir;

# ------------------------------------------------------------
# Config parameters
# ------------------------------------------------------------
my %CONF = read_config($config);

if (exists $CONF{ANALYZE_CORE_GENES_ONLY}) {
    die "[ERROR] ANALYZE_CORE_GENES_ONLY has been removed.\n"
      . "        run_gene_mean_depth.pl now always analyzes HCM core genes from HCM_CORE_GENE_LIST.\n"
      . "        Please delete ANALYZE_CORE_GENES_ONLY from the config file.\n";
}

my $project_root = detect_project_root($config);

my $samtools = get_conf_required(\%CONF, "SAMTOOLS");
check_executable($samtools, "SAMTOOLS");

my $gene_txt = resolve_config_path(
    get_conf_required(\%CONF, "REFSEQ_MANE_SELECT_GENE_TXT"),
    $project_root
);

die "[ERROR] REFSEQ_MANE_SELECT_GENE_TXT file not found: $gene_txt\n"
    unless -s $gene_txt;

my $core_gene_list = resolve_config_path(
    get_conf_required(\%CONF, "HCM_CORE_GENE_LIST"),
    $project_root
);

die "[ERROR] HCM_CORE_GENE_LIST file not found: $core_gene_list\n"
    unless -s $core_gene_list;

my $window_size = get_conf_required(\%CONF, "WINDOW_SIZE");
my $window_step = get_conf_required(\%CONF, "WINDOW_STEP");
my $min_window_size = get_conf_required(\%CONF, "MIN_WINDOW_SIZE");
my $min_gene_mean_depth = get_conf_required(\%CONF, "MIN_GENE_MEAN_DEPTH");
my $min_mapq = get_conf_required(\%CONF, "MIN_MAPQ");
my $exclude_dup = get_conf_required(\%CONF, "EXCLUDE_DUPLICATES");
my $del_depth_ratio_cutoff = get_conf_required(\%CONF, "DEL_DEPTH_RATIO_CUTOFF");
my $min_consecutive_del_windows = get_conf_required(\%CONF, "MIN_CONSECUTIVE_DEL_WINDOWS");
my $keep_depth_file = get_conf_required(\%CONF, "KEEP_GENE_DEPTH_FILE");

validate_positive_integer("WINDOW_SIZE", $window_size);
validate_positive_integer("WINDOW_STEP", $window_step);
validate_positive_integer("MIN_WINDOW_SIZE", $min_window_size);
validate_nonnegative_number("MIN_GENE_MEAN_DEPTH", $min_gene_mean_depth);
validate_nonnegative_integer("MIN_MAPQ", $min_mapq);
validate_bool_01("EXCLUDE_DUPLICATES", $exclude_dup);
validate_fraction_0_1("DEL_DEPTH_RATIO_CUTOFF", $del_depth_ratio_cutoff);
validate_positive_integer("MIN_CONSECUTIVE_DEL_WINDOWS", $min_consecutive_del_windows);
validate_bool_01("KEEP_GENE_DEPTH_FILE", $keep_depth_file);

# ------------------------------------------------------------
# Output files
# ------------------------------------------------------------
my $prefix = $out;
$prefix =~ s/\.tsv$//i;

my $depth_dir = $prefix . ".gene_depth_files";
make_path($depth_dir) unless -d $depth_dir;

my $gene_depth_out  = $prefix . ".gene_mean_depth.tsv";
my $window_depth_out = $prefix . ".window_depth.tsv";
my $all_ratio_out   = $prefix . ".all_window_ratio.tsv";
my $del_window_out  = $prefix . ".del_windows.tsv";

# ------------------------------------------------------------
# Read core gene list and gene regions
# ------------------------------------------------------------
my ($core_gene_ref, $classification_ref) = read_core_gene_list($core_gene_list);
my @genes = read_core_gene_regions(
    gene_txt           => $gene_txt,
    core_gene_ref      => $core_gene_ref,
    classification_ref => $classification_ref,
);

die "[ERROR] No core gene regions found in REFSEQ_MANE_SELECT_GENE_TXT: $gene_txt\n"
    unless @genes;

my %observed_gene = map { $_->{gene} => 1 } @genes;
my @missing_genes = grep { !exists $observed_gene{$_} } sort keys %{$core_gene_ref};

if (@missing_genes) {
    die "[ERROR] The following genes in HCM_CORE_GENE_LIST were not found in REFSEQ_MANE_SELECT_GENE_TXT:\n"
      . join("\n", map { "        $_" } @missing_genes)
      . "\n";
}

print "[INFO] Window-based depth analysis started\n";
print "[INFO] Sample : $sample\n";
print "[INFO] BAM : $bam\n";
print "[INFO] Config : $config\n";
print "[INFO] Gene coordinate file : $gene_txt\n";
print "[INFO] HCM core gene list : $core_gene_list\n";
print "[INFO] HCM core genes loaded : " . scalar(keys %{$core_gene_ref}) . "\n";
print "[INFO] Gene regions retained : " . scalar(@genes) . "\n";
print "[INFO] WINDOW_SIZE : $window_size\n";
print "[INFO] WINDOW_STEP : $window_step\n";
print "[INFO] MIN_WINDOW_SIZE : $min_window_size\n";
print "[INFO] MIN_GENE_MEAN_DEPTH : $min_gene_mean_depth\n";
print "[INFO] MIN_MAPQ : $min_mapq\n";
print "[INFO] EXCLUDE_DUPLICATES : $exclude_dup\n";
print "[INFO] DEL_DEPTH_RATIO_CUTOFF : $del_depth_ratio_cutoff\n";
print "[INFO] MIN_CONSECUTIVE_DEL_WINDOWS : $min_consecutive_del_windows\n";
print "[INFO] KEEP_GENE_DEPTH_FILE : $keep_depth_file\n";

# ------------------------------------------------------------
# Generate windows
# ------------------------------------------------------------
my @windows;
my %gene_windows;

for my $gene (@genes) {
    my $gkey = gene_key($gene);

    my @gene_windows = generate_windows(
        gene            => $gene,
        sample          => $sample,
        window_size     => $window_size,
        window_step     => $window_step,
        min_window_size => $min_window_size,
    );

    if (!@gene_windows) {
        warn "[WARNING] No window generated for gene "
          . $gene->{gene} . " "
          . $gene->{transcript} . " "
          . $gene->{chr} . ":" . $gene->{start} . "-" . $gene->{end}
          . ". Please check WINDOW_SIZE/MIN_WINDOW_SIZE.\n";
    }

    push @windows, @gene_windows;
    $gene_windows{$gkey} = \@gene_windows;
}

die "[ERROR] No windows generated. Please check gene coordinates and window parameters.\n"
    unless @windows;

print "[INFO] Total windows generated : " . scalar(@windows) . "\n";

# ------------------------------------------------------------
# Generate or reuse per-gene depth files and calculate depths
# ------------------------------------------------------------
my %gene_mean_depth;
my %gene_total_length;
my %window_mean_depth;

my $reused_depth_file_count = 0;
my $generated_depth_file_count = 0;

my $gene_index = 0;
my $gene_total = scalar(@genes);

for my $gene (@genes) {
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
      . $gene->{chr} . ":" . $gene->{start} . "-" . $gene->{end}
      . "\n";

    if (-s $depth_file) {
        print "[INFO] [$gene_index/$gene_total] Reuse existing gene depth file: $depth_file\n";
        $reused_depth_file_count++;
    } else {
        print "[INFO] [$gene_index/$gene_total] Generate gene depth file: $depth_file\n";
        run_samtools_depth_for_gene(
            samtools    => $samtools,
            bam         => $bam,
            gene        => $gene,
            depth_file  => $depth_file,
            min_mapq    => $min_mapq,
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

    $gene_mean_depth{$gkey} = $gene_mean;
    $gene_total_length{$gkey} = $gene_len;

    for my $wkey (keys %{$win_depth_ref}) {
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

# ------------------------------------------------------------
# Write gene mean depth output
# ------------------------------------------------------------
open my $GD, ">", $gene_depth_out
    or die "[ERROR] Cannot write gene mean depth output: $gene_depth_out\n";

print $GD join("\t", qw(
    SampleID
    Gene
    Classification
    Transcript
    Chrom
    Start
    End
    Strand
    ExonCount
    Gene_Key
    Gene_Length
    Gene_Mean_Depth
    Window_Count
)), "\n";

for my $gene (@genes) {
    my $gkey = gene_key($gene);
    my $window_count = exists $gene_windows{$gkey} ? scalar(@{ $gene_windows{$gkey} }) : 0;

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

# ------------------------------------------------------------
# Write window depth output
# ------------------------------------------------------------
open my $WD, ">", $window_depth_out
    or die "[ERROR] Cannot write window depth output: $window_depth_out\n";

print $WD join("\t", qw(
    SampleID
    Gene
    Classification
    Transcript
    Window_ID
    Chrom
    Start
    End
    Length
    Gene_Start
    Gene_End
    Window_Mean_Depth
)), "\n";

for my $win (@windows) {
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

# ------------------------------------------------------------
# Calculate depth ratios and classify windows
# ------------------------------------------------------------
my @ratio_records;
my $del_window_count = 0;

open my $RATIO, ">", $all_ratio_out
    or die "[ERROR] Cannot write all window ratio output: $all_ratio_out\n";

open my $DELWIN, ">", $del_window_out
    or die "[ERROR] Cannot write deletion window output: $del_window_out\n";

my @ratio_header = qw(
    SampleID
    Gene
    Classification
    Transcript
    Window_ID
    Chrom
    Start
    End
    Length
    Gene_Start
    Gene_End
    Gene_Mean_Depth
    Window_Mean_Depth
    Depth_Ratio
    Is_Del_Window
    Window_Status
    Comment
);

print $RATIO join("\t", @ratio_header), "\n";
print $DELWIN join("\t", @ratio_header), "\n";

for my $win (@windows) {
    my $wkey = window_key($win);
    my $gkey = $win->{gene_key};

    my $gene_mean = $gene_mean_depth{$gkey};
    my $window_mean = $window_mean_depth{$wkey} // 0;

    my $ratio = "NA";
    $ratio = $window_mean / $gene_mean
        if defined $gene_mean && $gene_mean > 0;

    my ($status, $comment, $is_del_window) = classify_window_depth(
        gene_mean_depth    => $gene_mean,
        ratio              => $ratio,
        min_gene_mean_depth => $min_gene_mean_depth,
        del_ratio_cutoff   => $del_depth_ratio_cutoff,
    );

    my $record = {
        %{$win},
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
        defined $gene_mean ? sprintf("%.4f", $gene_mean) : "NA",
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

# ------------------------------------------------------------
# Merge consecutive deletion-supporting windows into candidates
# ------------------------------------------------------------
my @candidates = merge_consecutive_del_windows(
    records                         => \@ratio_records,
    min_consecutive_del_windows      => $min_consecutive_del_windows,
);

# ------------------------------------------------------------
# Write candidate output
# ------------------------------------------------------------
open my $CAND, ">", $out
    or die "[ERROR] Cannot write candidate output: $out\n";

print $CAND join("\t", qw(
    SampleID
    Gene
    Classification
    Transcript
    Chrom
    Start
    End
    Candidate_Length
    Window_Count
    Gene_Mean_Depth
    Candidate_Mean_Depth
    Mean_Depth_Ratio
    Min_Depth_Ratio
    Max_Depth_Ratio
    Candidate_Status
    Del_Window_IDs
    Comment
)), "\n";

for my $cand (@candidates) {
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
print "[INFO] Sample : $sample\n";
print "[INFO] Gene mean depth output : $gene_depth_out\n";
print "[INFO] Window depth output : $window_depth_out\n";
print "[INFO] All window ratio output : $all_ratio_out\n";
print "[INFO] Del window output : $del_window_out\n";
print "[INFO] Candidate output : $out\n";
print "[INFO] Gene depth files reused : $reused_depth_file_count\n";
print "[INFO] Gene depth files generated : $generated_depth_file_count\n";
print "[INFO] Del window number : $del_window_count\n";
print "[INFO] Candidate number : " . scalar(@candidates) . "\n";

exit 0;

# ============================================================
# Config helpers
# ============================================================
sub read_config {
    my ($file) = @_;

    my %conf;

    open my $FH, "<", $file
        or die "[ERROR] Cannot open config file: $file\n";

    while (my $line = <$FH>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        $line =~ s/\s+#.*$//;

        next unless $line =~ /^\s*([^=\s]+)\s*=\s*(.*?)\s*$/;

        my ($key, $val) = ($1, $2);

        $key =~ s/^\s+|\s+$//g;
        $val =~ s/^\s+|\s+$//g;
        $val =~ s/^["']//;
        $val =~ s/["']$//;

        die "[ERROR] Empty config key found in $file\n"
            if $key eq "";

        $conf{$key} = $val;
    }

    close $FH;

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

sub detect_project_root {
    my ($config_path) = @_;

    my $config_dir = dirname($config_path);

    if ($config_dir =~ m{/conf$}) {
        return abs_path(dirname($config_dir));
    }

    return abs_path($config_dir);
}

sub resolve_config_path {
    my ($path, $project_root) = @_;

    die "[ERROR] Empty config path provided\n"
        unless defined $path && $path ne "";

    if ($path =~ m{^/}) {
        my $abs = abs_path($path);
        die "[ERROR] Path not found: $path\n"
            unless defined $abs;
        return $abs;
    }

    my $full = "$project_root/$path";
    my $abs = abs_path($full);

    die "[ERROR] Path not found: $full\n"
        unless defined $abs;

    return $abs;
}

# ============================================================
# Validation helpers
# ============================================================
sub validate_bool_01 {
    my ($name, $v) = @_;

    die "[ERROR] $name must be 0 or 1. Observed: $v\n"
        unless defined $v && $v =~ /^[01]$/;

    return 1;
}

sub validate_positive_integer {
    my ($name, $v) = @_;

    die "[ERROR] $name must be a positive integer. Observed: $v\n"
        unless defined $v && $v =~ /^\d+$/ && $v >= 1;

    return 1;
}

sub validate_nonnegative_integer {
    my ($name, $v) = @_;

    die "[ERROR] $name must be a non-negative integer. Observed: $v\n"
        unless defined $v && $v =~ /^\d+$/ && $v >= 0;

    return 1;
}

sub validate_nonnegative_number {
    my ($name, $v) = @_;

    die "[ERROR] $name must be numeric. Observed: $v\n"
        unless defined $v && $v =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;

    die "[ERROR] $name must be >= 0. Observed: $v\n"
        unless $v >= 0;

    return 1;
}

sub validate_fraction_0_1 {
    my ($name, $v) = @_;

    die "[ERROR] $name must be numeric. Observed: $v\n"
        unless defined $v && $v =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;

    die "[ERROR] $name must be > 0 and <= 1. Observed: $v\n"
        unless $v > 0 && $v <= 1;

    return 1;
}

sub check_executable {
    my ($cmd, $name) = @_;

    die "[ERROR] Empty executable path for $name\n"
        unless defined $cmd && $cmd ne "";

    if ($cmd =~ m{/}) {
        die "[ERROR] $name executable not found or not executable: $cmd\n"
            unless -x $cmd;
    } else {
        my $ret = system("command -v $cmd >/dev/null 2>&1");
        die "[ERROR] $name executable not found in PATH: $cmd\n"
            if $ret != 0;
    }

    return 1;
}

sub check_bam_index {
    my ($bam_file) = @_;

    my $bai1 = "$bam_file.bai";
    my $bai2 = $bam_file;
    $bai2 =~ s/\.bam$/.bai/i;

    die "[ERROR] BAM index not found for: $bam_file\n"
      . "        Expected one of:\n"
      . "        $bai1\n"
      . "        $bai2\n"
        unless -s $bai1 || -s $bai2;

    return 1;
}

# ============================================================
# Input parsers
# ============================================================
sub read_core_gene_list {
    my ($file) = @_;

    my %gene;
    my %classification;

    open my $FH, "<", $file
        or die "[ERROR] Cannot open HCM_CORE_GENE_LIST: $file\n";

    my $line_no = 0;

    while (my $line = <$FH>) {
        chomp $line;
        $line =~ s/\r$//;
        $line =~ s/^\s+|\s+$//g;
        $line_no++;

        next if $line eq "";
        next if $line =~ /^#/;

        my @f = split /\t|\s+/, $line;

        my $g = $f[0] // "";
        my $c = $f[1] // "NA";

        $g =~ s/^\s+|\s+$//g;
        $c =~ s/^\s+|\s+$//g;

        next if $line_no == 1 && $g =~ /^Gene$/i;

        die "[ERROR] Empty gene symbol in HCM_CORE_GENE_LIST at line $line_no\n"
            if $g eq "";

        die "[ERROR] Invalid gene symbol in HCM_CORE_GENE_LIST at line $line_no: $g\n"
            unless $g =~ /^[A-Za-z0-9_.-]+$/;

        $c = "NA" if $c eq "";

        $gene{$g} = 1;

        if (!exists $classification{$g}) {
            $classification{$g} = $c;
        } else {
            $classification{$g} = choose_higher_classification($classification{$g}, $c);
        }
    }

    close $FH;

    die "[ERROR] No genes loaded from HCM_CORE_GENE_LIST: $file\n"
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

sub read_core_gene_regions {
    my %args = @_;

    my $file = $args{gene_txt};
    my $core_gene_ref = $args{core_gene_ref};
    my $classification_ref = $args{classification_ref};

    my @genes;
    my %col;

    open my $FH, "<", $file
        or die "[ERROR] Cannot open REFSEQ_MANE_SELECT_GENE_TXT: $file\n";

    my $line_no = 0;
    my $header_seen = 0;

    while (my $line = <$FH>) {
        chomp $line;
        $line =~ s/\r$//;
        $line_no++;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        my @f = split /\t/, $line, -1;

        if (!$header_seen) {
            for my $i (0 .. $#f) {
                my $name = $f[$i];
                $name =~ s/^\s+|\s+$//g;
                $col{$name} = $i;
            }

            for my $required (qw/Gene Transcript Chrom Start End Strand ExonCount/) {
                die "[ERROR] Required column '$required' not found in REFSEQ_MANE_SELECT_GENE_TXT: $file\n"
                    unless exists $col{$required};
            }

            $header_seen = 1;
            next;
        }

        my $gene       = $f[$col{Gene}];
        my $transcript = $f[$col{Transcript}];
        my $chr        = $f[$col{Chrom}];
        my $start      = $f[$col{Start}];
        my $end        = $f[$col{End}];
        my $strand     = $f[$col{Strand}];
        my $exon_count = $f[$col{ExonCount}];

        for ($gene, $transcript, $chr, $start, $end, $strand, $exon_count) {
            $_ = "" unless defined $_;
            s/^\s+|\s+$//g;
        }

        next unless exists $core_gene_ref->{$gene};

        die "[ERROR] Invalid Start at REFSEQ_MANE_SELECT_GENE_TXT line $line_no: $start\n"
            unless $start =~ /^\d+$/ && $start >= 1;

        die "[ERROR] Invalid End at REFSEQ_MANE_SELECT_GENE_TXT line $line_no: $end\n"
            unless $end =~ /^\d+$/ && $end >= 1;

        die "[ERROR] Start > End at REFSEQ_MANE_SELECT_GENE_TXT line $line_no: $chr:$start-$end\n"
            if $start > $end;

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

    die "[ERROR] REFSEQ_MANE_SELECT_GENE_TXT header was not found: $file\n"
        unless $header_seen;

    @genes = sort {
        chrom_order_key($a->{chr}) cmp chrom_order_key($b->{chr})
            || $a->{start} <=> $b->{start}
            || $a->{end} <=> $b->{end}
            || $a->{gene} cmp $b->{gene}
    } @genes;

    return @genes;
}

# ============================================================
# Window and depth calculation
# ============================================================
sub generate_windows {
    my %args = @_;

    my $gene = $args{gene};
    my $sample = $args{sample};
    my $window_size = $args{window_size};
    my $window_step = $args{window_step};
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

    my $samtools = $args{samtools};
    my $bam = $args{bam};
    my $gene = $args{gene};
    my $depth_file = $args{depth_file};
    my $min_mapq = $args{min_mapq};
    my $exclude_dup = $args{exclude_dup};

    my $region = join_region($gene->{chr}, $gene->{start}, $gene->{end});

    # Explicitly discard low-confidence alignment categories.
    # 0x4   read unmapped
    # 0x100 secondary alignment
    # 0x200 QC fail
    # 0x800 supplementary alignment
    # 0x400 duplicate, only when EXCLUDE_DUPLICATES=1
    my $exclude_flags = 4 + 256 + 512 + 2048;
    $exclude_flags += 1024 if $exclude_dup == 1;

    my @cmd = (
        $samtools, "depth",
        "-a",
        "-Q", $min_mapq,
        "-G", $exclude_flags,
        "-r", $region,
        $bam,
    );

    print "[INFO] Generate gene depth file: $depth_file\n";
    print "[CMD] ", join(" ", map { shell_quote($_) } @cmd), " > ", shell_quote($depth_file), "\n";

    open my $IN, "-|", @cmd
        or die "[ERROR] Failed to run samtools depth for region $region\n";

    open my $OUT, ">", $depth_file
        or die "[ERROR] Cannot write depth file: $depth_file\n";

    while (my $line = <$IN>) {
        print $OUT $line;
    }

    close $IN
        or die "[ERROR] samtools depth failed for region $region\n";

    close $OUT;

    die "[ERROR] Depth file was not generated: $depth_file\n"
        unless -e $depth_file;
}

sub calculate_depth_from_gene_depth_file {
    my %args = @_;

    my $depth_file = $args{depth_file};
    my $gene = $args{gene};
    my $windows_ref = $args{windows_ref};
    my @windows = @{$windows_ref};

    my $gene_start = $gene->{start};
    my $gene_end = $gene->{end};
    my $gene_length = $gene_end - $gene_start + 1;

    die "[ERROR] Invalid gene length for " . $gene->{gene} . ": $gene_start-$gene_end\n"
        if $gene_length <= 0;

    my @depth_by_offset;

    open my $IN, "<", $depth_file
        or die "[ERROR] Cannot open depth file: $depth_file\n";

    while (my $line = <$IN>) {
        chomp $line;
        next if $line =~ /^\s*$/;

        my @f = split /\t/, $line;
        next unless @f >= 3;

        my ($chr, $pos, $depth) = @f[0, 1, 2];

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

    my $gene_mean_depth = $prefix_sum[$gene_length] / $gene_length;

    my %window_mean_depth;

    for my $win (@windows) {
        my $wkey = window_key($win);

        my $w_start = $win->{start};
        my $w_end = $win->{end};

        $w_start = $gene_start if $w_start < $gene_start;
        $w_end = $gene_end if $w_end > $gene_end;

        my $w_len = $w_end - $w_start + 1;

        if ($w_len <= 0) {
            $window_mean_depth{$wkey} = 0;
            next;
        }

        my $left_offset = $w_start - $gene_start + 1;
        my $right_offset = $w_end - $gene_start + 1;

        my $window_depth_sum = $prefix_sum[$right_offset] - $prefix_sum[$left_offset - 1];
        my $mean = $window_depth_sum / $w_len;

        $window_mean_depth{$wkey} = $mean;
    }

    return ($gene_mean_depth, $gene_length, \%window_mean_depth);
}

sub classify_window_depth {
    my %args = @_;

    my $gene_mean_depth = $args{gene_mean_depth};
    my $ratio = $args{ratio};
    my $min_gene_mean_depth = $args{min_gene_mean_depth};
    my $cutoff = $args{del_ratio_cutoff};

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

# ============================================================
# Candidate merging
# ============================================================
sub merge_consecutive_del_windows {
    my %args = @_;

    my @records = @{ $args{records} };
    my $min_consecutive = $args{min_consecutive_del_windows};

    my %grouped;

    for my $r (@records) {
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

    for my $gkey (sort keys %grouped) {
        my @sorted = sort {
            $a->{start} <=> $b->{start}
                || $a->{end} <=> $b->{end}
        } @{ $grouped{$gkey} };

        my @current_del_run;

        for my $r (@sorted) {
            if ($r->{is_del_window}) {
                push @current_del_run, $r;
            } else {
                add_del_run_candidate_if_valid(
                    candidates_ref => \@candidates,
                    windows_ref    => \@current_del_run,
                    min_consecutive => $min_consecutive,
                );
                @current_del_run = ();
            }
        }

        add_del_run_candidate_if_valid(
            candidates_ref => \@candidates,
            windows_ref    => \@current_del_run,
            min_consecutive => $min_consecutive,
        );
    }

    return @candidates;
}

sub add_del_run_candidate_if_valid {
    my %args = @_;

    my $candidates_ref = $args{candidates_ref};
    my $windows_ref = $args{windows_ref};
    my $min_consecutive = $args{min_consecutive};

    return unless @{$windows_ref};
    return if scalar(@{$windows_ref}) < $min_consecutive;

    my $candidate = build_candidate($windows_ref);
    push @{$candidates_ref}, $candidate;
}

sub build_candidate {
    my ($wins_ref) = @_;

    my @wins = sort {
        $a->{start} <=> $b->{start}
            || $a->{end} <=> $b->{end}
    } @{$wins_ref};

    my $first = $wins[0];
    my $last = $wins[-1];

    my $start = $first->{start};
    my $end = $last->{end};

    my $total_depth_sum = 0;
    my $total_length = 0;
    my $ratio_sum = 0;
    my $ratio_count = 0;
    my $min_ratio;
    my $max_ratio;
    my @window_ids;

    for my $w (@wins) {
        push @window_ids, $w->{window_id};

        $total_depth_sum += $w->{window_mean_depth} * $w->{length};
        $total_length += $w->{length};

        if (defined $w->{depth_ratio} && $w->{depth_ratio} ne "NA") {
            $ratio_sum += $w->{depth_ratio};
            $ratio_count++;

            $min_ratio = $w->{depth_ratio}
                if !defined $min_ratio || $w->{depth_ratio} < $min_ratio;

            $max_ratio = $w->{depth_ratio}
                if !defined $max_ratio || $w->{depth_ratio} > $max_ratio;
        }
    }

    my $candidate_mean_depth = 0;
    $candidate_mean_depth = $total_depth_sum / $total_length
        if $total_length > 0;

    my $mean_depth_ratio = 0;
    $mean_depth_ratio = $ratio_sum / $ratio_count
        if $ratio_count > 0;

    $min_ratio = 0 unless defined $min_ratio;
    $max_ratio = 0 unless defined $max_ratio;

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
        min_depth_ratio      => $min_ratio,
        max_depth_ratio      => $max_ratio,
        del_window_ids       => join(",", @window_ids),
        comment              => "Merged consecutive deletion-supporting windows",
    };
}

# ============================================================
# Utility helpers
# ============================================================
sub make_gene_depth_filename {
    my %args = @_;

    my $depth_dir = $args{depth_dir};
    my $sample = $args{sample};
    my $gene = $args{gene};

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
    my $end = exists $r->{gene_end} ? $r->{gene_end} : $r->{end};

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

sub chrom_order_key {
    my ($chr) = @_;

    my $x = $chr;
    $x =~ s/^chr//i;

    return sprintf("%03d", $x) if $x =~ /^\d+$/;

    my $u = uc($x);

    return "023" if $u eq "X";
    return "024" if $u eq "Y";
    return "025" if $u eq "M" || $u eq "MT";

    return "999_$chr";
}

sub shell_quote {
    my ($str) = @_;

    die "[ERROR] Undefined shell argument\n"
        unless defined $str;

    return "''" if $str eq "";

    return $str
        if $str =~ /^[A-Za-z0-9_\.\-\/\:=,\+]+$/;

    $str =~ s/'/'\\''/g;

    return "'$str'";
}

sub usage {
    return <<'USAGE';
Usage:
  perl bin/run_gene_mean_depth.pl \
    --config /abs/path/conf/hcm_exondel.conf \
    --bam /abs/path/SAMPLE.target.bam \
    --sample SAMPLE_ID \
    --out /abs/path/results/SAMPLE/01.depth/SAMPLE.depth_candidates.tsv

Required arguments:
  --config    HCMExonDel config file.
  --bam       Coordinate-sorted and indexed target BAM generated by extract_target_bam.pl.
  --sample    Sample ID.
  --out       Output depth candidate TSV file.

Required config parameters:
  SAMTOOLS
  REFSEQ_MANE_SELECT_GENE_TXT
  HCM_CORE_GENE_LIST
  WINDOW_SIZE
  WINDOW_STEP
  MIN_WINDOW_SIZE
  MIN_GENE_MEAN_DEPTH
  MIN_MAPQ
  EXCLUDE_DUPLICATES
  DEL_DEPTH_RATIO_CUTOFF
  MIN_CONSECUTIVE_DEL_WINDOWS
  KEEP_GENE_DEPTH_FILE

Outputs:
  Main output:
    *.depth_candidates.tsv

  Additional outputs:
    *.gene_mean_depth.tsv
    *.window_depth.tsv
    *.all_window_ratio.tsv
    *.del_windows.tsv

  Intermediate output:
    *.gene_depth_files/*.depth.tsv

Notes:
  1. HCMExonDel only analyzes genes listed in HCM_CORE_GENE_LIST.
  2. This script is intended to run on SAMPLE.target.bam, not the original WGS BAM.
  3. If MIN_MAPQ, EXCLUDE_DUPLICATES, BAM, or gene coordinates are changed,
     delete the old *.gene_depth_files directory before rerunning.
USAGE
}

