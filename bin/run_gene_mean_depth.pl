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
# Window-based depth ratio analysis for HCM WGS exon deletion
# detection.
#
# Main idea:
#   1. Split target regions into fixed-size windows
#   2. Calculate mean depth for each window
#   3. Calculate region/gene mean depth
#   4. Calculate Window_Mean_Depth / Region_Mean_Depth
#   5. Merge consecutive low-depth windows as deletion candidates
#
# Expected target BED format:
#
#   chr18    33877675    34360183    FHOD3|NM_001281740.3|gene
#
# or:
#
#   chr18    34232239    34241309    FHOD3|NM_001281740.3|EX12_EX14
#
# BED is 0-based start, 1-based end.
# Output is 1-based start and 1-based end.
# ============================================================

my $config;
my $bam;
my $sample;
my $target_bed;
my $out;
my $threads = 4;
my $window_size;
my $window_step;
my $help = 0;

GetOptions(
    "config=s"     => \$config,
    "bam=s"        => \$bam,
    "sample=s"     => \$sample,
    "target-bed=s" => \$target_bed,
    "out=s"        => \$out,
    "threads=i"    => \$threads,
    "window-size=i" => \$window_size,
    "window-step=i" => \$window_step,
    "help"         => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless $config && $bam && $sample && $target_bed && $out;

$config     = abs_path($config);
$bam        = abs_path($bam);
$target_bed = abs_path($target_bed);

die "[ERROR] Config file not found: $config\n" unless -s $config;
die "[ERROR] BAM file not found: $bam\n" unless -s $bam;
die "[ERROR] Target BED file not found: $target_bed\n" unless -s $target_bed;

my %CONF = read_config($config);

my $samtools = $CONF{SAMTOOLS} || "samtools";

my $min_mapq = defined $CONF{MIN_MAPQ} ? $CONF{MIN_MAPQ} : 20;

$window_size ||= $CONF{WINDOW_SIZE} || 1000;
$window_step ||= $CONF{WINDOW_STEP} || $window_size;

my $min_region_mean_depth = defined $CONF{MIN_REGION_MEAN_DEPTH}
    ? $CONF{MIN_REGION_MEAN_DEPTH}
    : defined $CONF{MIN_GENE_MEAN_DEPTH}
        ? $CONF{MIN_GENE_MEAN_DEPTH}
        : 20;

my $depth_ratio_cutoff = defined $CONF{DEPTH_RATIO_CUTOFF}
    ? $CONF{DEPTH_RATIO_CUTOFF}
    : 0.65;

my $het_del_ratio_low = defined $CONF{HET_DEL_RATIO_LOW}
    ? $CONF{HET_DEL_RATIO_LOW}
    : 0.35;

my $het_del_ratio_high = defined $CONF{HET_DEL_RATIO_HIGH}
    ? $CONF{HET_DEL_RATIO_HIGH}
    : 0.70;

my $min_candidate_windows = defined $CONF{MIN_CANDIDATE_WINDOWS}
    ? $CONF{MIN_CANDIDATE_WINDOWS}
    : 1;

my $max_merge_gap_windows = defined $CONF{MAX_MERGE_GAP_WINDOWS}
    ? $CONF{MAX_MERGE_GAP_WINDOWS}
    : 0;

die "[ERROR] WINDOW_SIZE must be a positive integer\n" unless $window_size > 0;
die "[ERROR] WINDOW_STEP must be a positive integer\n" unless $window_step > 0;

my $outdir = dirname($out);
make_path($outdir) unless -d $outdir;

my $prefix = $out;
$prefix =~ s/\.tsv$//;

my $window_depth_out = $prefix . ".window_depth.tsv";
my $region_depth_out = $prefix . ".region_mean_depth.tsv";
my $all_ratio_out    = $prefix . ".all_window_ratio.tsv";

# ============================================================
# Step 1. Read target regions
# ============================================================

my @regions;

open my $bed_fh, "<", $target_bed
    or die "[ERROR] Cannot open target BED: $target_bed\n";

while (my $line = <$bed_fh>) {
    chomp $line;
    $line =~ s/\r$//;

    next if $line =~ /^\s*$/;
    next if $line =~ /^\s*#/;

    my @f = split /\t/, $line;
    die "[ERROR] Invalid BED line: $line\n" unless @f >= 4;

    my ($chr, $start0, $end1, $name) = @f[0, 1, 2, 3];

    die "[ERROR] Invalid BED coordinate: $line\n"
        unless $start0 =~ /^\d+$/ && $end1 =~ /^\d+$/ && $end1 > $start0;

    my ($gene, $transcript, $region_name) = parse_region_name($name);

    push @regions, {
        chr         => $chr,
        start0      => $start0,
        start1      => $start0 + 1,
        end1        => $end1,
        length      => $end1 - $start0,
        name        => $name,
        gene        => $gene,
        transcript  => $transcript,
        region_name => $region_name,
    };
}

close $bed_fh;

die "[ERROR] No valid target regions found in BED: $target_bed\n" unless @regions;

# ============================================================
# Step 2. Generate fixed windows
# ============================================================

my @windows;
my %region_windows;

foreach my $region (@regions) {
    my $region_key = region_key($region);
    my $win_index  = 0;

    for (my $s0 = $region->{start0}; $s0 < $region->{end1}; $s0 += $window_step) {
        my $e1 = $s0 + $window_size;
        $e1 = $region->{end1} if $e1 > $region->{end1};

        next if $e1 <= $s0;

        $win_index++;

        my $window = {
            sample      => $sample,
            chr         => $region->{chr},
            start0      => $s0,
            start1      => $s0 + 1,
            end1        => $e1,
            length      => $e1 - $s0,
            gene        => $region->{gene},
            transcript  => $region->{transcript},
            region_name => $region->{region_name},
            region_key  => $region_key,
            window_id   => $region->{gene} . "_W" . $win_index,
        };

        push @windows, $window;
        push @{ $region_windows{$region_key} }, $window;
    }
}

die "[ERROR] No windows generated. Please check target BED and WINDOW_SIZE\n" unless @windows;

# ============================================================
# Step 3. Calculate mean depth for each window
# ============================================================

my %window_depth;

foreach my $win (@windows) {
    my $region = $win->{chr} . ":" . $win->{start1} . "-" . $win->{end1};

    my $mean_depth = calculate_region_mean_depth(
        samtools => $samtools,
        bam      => $bam,
        region   => $region,
        length   => $win->{length},
        min_mapq => $min_mapq,
    );

    my $key = window_key($win);
    $window_depth{$key} = $mean_depth;
}

# ============================================================
# Step 4. Calculate region mean depth
# ============================================================

my %region_mean_depth;
my %region_total_length;

foreach my $region_key (sort keys %region_windows) {
    my $total_depth_sum = 0;
    my $total_length    = 0;

    foreach my $win (@{ $region_windows{$region_key} }) {
        my $key = window_key($win);
        my $mean_depth = $window_depth{$key};

        $total_depth_sum += $mean_depth * $win->{length};
        $total_length    += $win->{length};
    }

    my $region_mean = 0;
    if ($total_length > 0) {
        $region_mean = $total_depth_sum / $total_length;
    }

    $region_mean_depth{$region_key} = $region_mean;
    $region_total_length{$region_key} = $total_length;
}

# ============================================================
# Step 5. Output window depth and region mean depth
# ============================================================

open my $wd_fh, ">", $window_depth_out
    or die "[ERROR] Cannot write window depth output: $window_depth_out\n";

print $wd_fh join("\t",
    qw(
        SampleID
        Gene
        Transcript
        Region
        Window_ID
        Chrom
        Start
        End
        Length
        Window_Mean_Depth
    )
) . "\n";

foreach my $win (@windows) {
    my $key = window_key($win);

    print $wd_fh join("\t",
        $sample,
        $win->{gene},
        $win->{transcript},
        $win->{region_name},
        $win->{window_id},
        $win->{chr},
        $win->{start1},
        $win->{end1},
        $win->{length},
        sprintf("%.4f", $window_depth{$key}),
    ) . "\n";
}

close $wd_fh;


open my $rd_fh, ">", $region_depth_out
    or die "[ERROR] Cannot write region mean depth output: $region_depth_out\n";

print $rd_fh join("\t",
    qw(
        SampleID
        Gene
        Transcript
        Region
        Region_Key
        Region_Mean_Depth
        Window_Count
        Total_Window_Length
    )
) . "\n";

foreach my $region_key (sort keys %region_windows) {
    my $first = $region_windows{$region_key}->[0];

    print $rd_fh join("\t",
        $sample,
        $first->{gene},
        $first->{transcript},
        $first->{region_name},
        $region_key,
        sprintf("%.4f", $region_mean_depth{$region_key}),
        scalar(@{ $region_windows{$region_key} }),
        $region_total_length{$region_key},
    ) . "\n";
}

close $rd_fh;

# ============================================================
# Step 6. Calculate window ratio and mark low-depth windows
# ============================================================

my @ratio_records;

open my $ratio_fh, ">", $all_ratio_out
    or die "[ERROR] Cannot write all window ratio output: $all_ratio_out\n";

print $ratio_fh join("\t",
    qw(
        SampleID
        Gene
        Transcript
        Region
        Window_ID
        Chrom
        Start
        End
        Length
        Region_Mean_Depth
        Window_Mean_Depth
        Depth_Ratio
        Window_Status
        Comment
    )
) . "\n";

foreach my $win (@windows) {
    my $wkey = window_key($win);
    my $rkey = $win->{region_key};

    my $region_mean = $region_mean_depth{$rkey};
    my $window_mean = $window_depth{$wkey};

    my $ratio = "NA";
    if ($region_mean && $region_mean > 0) {
        $ratio = $window_mean / $region_mean;
    }

    my ($status, $comment) = classify_window_depth(
        region_mean           => $region_mean,
        window_mean           => $window_mean,
        ratio                 => $ratio,
        min_region_mean_depth => $min_region_mean_depth,
        depth_ratio_cutoff    => $depth_ratio_cutoff,
        het_del_ratio_low     => $het_del_ratio_low,
        het_del_ratio_high    => $het_del_ratio_high,
    );

    my $record = {
        %$win,
        region_mean   => $region_mean,
        window_mean   => $window_mean,
        depth_ratio   => $ratio,
        window_status => $status,
        comment       => $comment,
    };

    push @ratio_records, $record;

    print $ratio_fh join("\t",
        $sample,
        $win->{gene},
        $win->{transcript},
        $win->{region_name},
        $win->{window_id},
        $win->{chr},
        $win->{start1},
        $win->{end1},
        $win->{length},
        sprintf("%.4f", $region_mean),
        sprintf("%.4f", $window_mean),
        ($ratio eq "NA" ? "NA" : sprintf("%.4f", $ratio)),
        $status,
        $comment,
    ) . "\n";
}

close $ratio_fh;

# ============================================================
# Step 7. Merge consecutive low-depth windows into candidates
# ============================================================

my @candidates = merge_low_depth_windows(
    records               => \@ratio_records,
    min_candidate_windows => $min_candidate_windows,
    max_merge_gap_windows => $max_merge_gap_windows,
);

open my $cand_fh, ">", $out
    or die "[ERROR] Cannot write candidate output: $out\n";

print $cand_fh join("\t",
    qw(
        SampleID
        Gene
        Transcript
        Region
        Chrom
        Start
        End
        Candidate_Length
        Window_Count
        Region_Mean_Depth
        Candidate_Mean_Depth
        Mean_Depth_Ratio
        Min_Depth_Ratio
        Max_Depth_Ratio
        Depth_Status
        Candidate_Status
        Comment
    )
) . "\n";

foreach my $cand (@candidates) {
    print $cand_fh join("\t",
        $sample,
        $cand->{gene},
        $cand->{transcript},
        $cand->{region_name},
        $cand->{chr},
        $cand->{start1},
        $cand->{end1},
        $cand->{candidate_length},
        $cand->{window_count},
        sprintf("%.4f", $cand->{region_mean}),
        sprintf("%.4f", $cand->{candidate_mean_depth}),
        sprintf("%.4f", $cand->{mean_depth_ratio}),
        sprintf("%.4f", $cand->{min_depth_ratio}),
        sprintf("%.4f", $cand->{max_depth_ratio}),
        $cand->{depth_status},
        "Depth_candidate",
        $cand->{comment},
    ) . "\n";
}

close $cand_fh;

print "[INFO] Window-based depth analysis finished\n";
print "[INFO] Window depth output       : $window_depth_out\n";
print "[INFO] Region mean depth output  : $region_depth_out\n";
print "[INFO] All window ratio output   : $all_ratio_out\n";
print "[INFO] Candidate output          : $out\n";
print "[INFO] Candidate number          : " . scalar(@candidates) . "\n";

exit 0;


# ============================================================
# Subroutines
# ============================================================

sub calculate_region_mean_depth {
    my %args = @_;

    my $samtools = $args{samtools};
    my $bam      = $args{bam};
    my $region   = $args{region};
    my $length   = $args{length};
    my $min_mapq = $args{min_mapq};

    my $cmd = join(" ",
        shell_quote($samtools),
        "depth",
        "-a",
        "-Q", $min_mapq,
        "-r", shell_quote($region),
        shell_quote($bam)
    );

    open my $pipe, "$cmd |"
        or die "[ERROR] Failed to run command: $cmd\n";

    my $depth_sum = 0;

    while (my $line = <$pipe>) {
        chomp $line;
        next if $line =~ /^\s*$/;

        my @f = split /\t/, $line;
        next unless @f >= 3;

        my $depth = $f[2];
        next unless defined $depth && $depth =~ /^\d+$/;

        $depth_sum += $depth;
    }

    close $pipe;

    my $mean_depth = 0;
    if ($length > 0) {
        $mean_depth = $depth_sum / $length;
    }

    return $mean_depth;
}


sub parse_region_name {
    my ($name) = @_;

    my ($gene, $transcript, $region_name);

    if ($name =~ /\|/) {
        my @x = split /\|/, $name;
        $gene        = $x[0] || "NA";
        $transcript  = $x[1] || "NA";
        $region_name = $x[2] || "target";
    }
    elsif ($name =~ /:/) {
        my @x = split /:/, $name;
        $gene        = $x[0] || "NA";
        $transcript  = $x[1] || "NA";
        $region_name = $x[2] || "target";
    }
    else {
        $gene        = $name;
        $transcript  = "NA";
        $region_name = "target";
    }

    return ($gene, $transcript, $region_name);
}


sub region_key {
    my ($region) = @_;

    return join("|",
        $region->{gene},
        $region->{transcript},
        $region->{region_name},
        $region->{chr},
        $region->{start1},
        $region->{end1},
    );
}


sub window_key {
    my ($win) = @_;

    return join("|",
        $win->{gene},
        $win->{transcript},
        $win->{region_name},
        $win->{window_id},
        $win->{chr},
        $win->{start1},
        $win->{end1},
    );
}


sub classify_window_depth {
    my %args = @_;

    my $region_mean = $args{region_mean};
    my $window_mean = $args{window_mean};
    my $ratio       = $args{ratio};

    my $min_region_mean_depth = $args{min_region_mean_depth};
    my $cutoff                = $args{depth_ratio_cutoff};
    my $het_low               = $args{het_del_ratio_low};
    my $het_high              = $args{het_del_ratio_high};

    if (!$region_mean || $region_mean < $min_region_mean_depth) {
        return (
            "Low_region_depth",
            "Region mean depth is lower than threshold"
        );
    }

    if ($ratio eq "NA") {
        return (
            "NA",
            "Depth ratio cannot be calculated"
        );
    }

    if ($ratio <= $cutoff) {
        if ($ratio >= $het_low && $ratio <= $het_high) {
            return (
                "Low_depth_window",
                "Depth ratio is consistent with possible heterozygous deletion"
            );
        }
        elsif ($ratio < $het_low) {
            return (
                "Very_low_depth_window",
                "Depth ratio is lower than expected heterozygous deletion range"
            );
        }
        else {
            return (
                "Mild_low_depth_window",
                "Depth ratio is below cutoff"
            );
        }
    }

    return (
        "Normal_depth_window",
        "Depth ratio is not reduced"
    );
}


sub merge_low_depth_windows {
    my %args = @_;

    my $records_ref = $args{records};
    my $min_windows = $args{min_candidate_windows};
    my $max_gap     = $args{max_merge_gap_windows};

    my @records = @$records_ref;

    my %grouped;

    foreach my $r (@records) {
        my $key = join("|",
            $r->{gene},
            $r->{transcript},
            $r->{region_name},
            $r->{chr},
            $r->{region_key},
        );

        push @{ $grouped{$key} }, $r;
    }

    my @candidates;

    foreach my $gkey (sort keys %grouped) {
        my @sorted = sort {
            $a->{start1} <=> $b->{start1}
        } @{ $grouped{$gkey} };

        my @current;
        my $gap_count = 0;

        for (my $i = 0; $i < @sorted; $i++) {
            my $r = $sorted[$i];

            my $is_low = is_low_depth_status($r->{window_status});

            if ($is_low) {
                push @current, $r;
                $gap_count = 0;
            }
            else {
                if (@current && $gap_count < $max_gap) {
                    $gap_count++;
                    push @current, $r;
                }
                else {
                    if (@current) {
                        my @candidate_windows = grep {
                            is_low_depth_status($_->{window_status})
                        } @current;

                        if (@candidate_windows >= $min_windows) {
                            push @candidates, build_candidate(\@current);
                        }

                        @current = ();
                        $gap_count = 0;
                    }
                }
            }
        }

        if (@current) {
            my @candidate_windows = grep {
                is_low_depth_status($_->{window_status})
            } @current;

            if (@candidate_windows >= $min_windows) {
                push @candidates, build_candidate(\@current);
            }
        }
    }

    return @candidates;
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

    my @wins = @$wins_ref;

    my @low_wins = grep {
        is_low_depth_status($_->{window_status})
    } @wins;

    my $first = $wins[0];
    my $last  = $wins[-1];

    my $start1 = $first->{start1};
    my $end1   = $last->{end1};

    my $total_depth_sum = 0;
    my $total_length    = 0;

    my $ratio_sum = 0;
    my $ratio_count = 0;

    my $min_ratio = 999999;
    my $max_ratio = -1;

    foreach my $w (@low_wins) {
        $total_depth_sum += $w->{window_mean} * $w->{length};
        $total_length    += $w->{length};

        if ($w->{depth_ratio} ne "NA") {
            $ratio_sum += $w->{depth_ratio};
            $ratio_count++;

            $min_ratio = $w->{depth_ratio} if $w->{depth_ratio} < $min_ratio;
            $max_ratio = $w->{depth_ratio} if $w->{depth_ratio} > $max_ratio;
        }
    }

    my $candidate_mean_depth = 0;
    if ($total_length > 0) {
        $candidate_mean_depth = $total_depth_sum / $total_length;
    }

    my $mean_depth_ratio = 0;
    if ($ratio_count > 0) {
        $mean_depth_ratio = $ratio_sum / $ratio_count;
    }

    my $depth_status = "Low_depth_region";
    my $comment = "Merged consecutive low-depth windows";

    if ($mean_depth_ratio >= 0.35 && $mean_depth_ratio <= 0.70) {
        $depth_status = "Heterozygous_deletion_like";
        $comment = "Merged low-depth windows are consistent with possible heterozygous deletion";
    }
    elsif ($mean_depth_ratio < 0.35) {
        $depth_status = "Severe_depth_reduction";
        $comment = "Merged low-depth windows show severe depth reduction";
    }

    return {
        gene                 => $first->{gene},
        transcript           => $first->{transcript},
        region_name          => $first->{region_name},
        chr                  => $first->{chr},
        start1               => $start1,
        end1                 => $end1,
        candidate_length     => $end1 - $start1 + 1,
        window_count         => scalar(@low_wins),
        region_mean          => $first->{region_mean},
        candidate_mean_depth => $candidate_mean_depth,
        mean_depth_ratio     => $mean_depth_ratio,
        min_depth_ratio      => $min_ratio == 999999 ? 0 : $min_ratio,
        max_depth_ratio      => $max_ratio == -1 ? 0 : $max_ratio,
        depth_status         => $depth_status,
        comment              => $comment,
    };
}


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


sub shell_quote {
    my ($str) = @_;

    return "''" unless defined $str;

    $str =~ s/'/'"'"'/g;
    return "'$str'";
}


sub usage {
    return <<"USAGE";

Usage:

  perl bin/run_gene_mean_depth.pl \\
    --config conf/hcm_exondel.conf \\
    --bam sample.sorted.bam \\
    --sample SAMPLE001 \\
    --target-bed db/target_region.bed \\
    --out results/SAMPLE001/01.depth/SAMPLE001.depth_candidates.tsv

Required arguments:

  --config       Config file
  --bam          Coordinate-sorted BAM file
  --sample       Sample ID
  --target-bed   Target region BED file
  --out          Output candidate file

Optional arguments:

  --threads       Number of threads [default: 4]
  --window-size   Window size in bp [default: config WINDOW_SIZE or 1000]
  --window-step   Window step in bp [default: config WINDOW_STEP or WINDOW_SIZE]
  --help          Show this help message

Expected target BED format:

  chr18    33877675    34360183    FHOD3|NM_001281740.3|gene

or:

  chr18    34232239    34241309    FHOD3|NM_001281740.3|EX12_EX14

Config parameters used by this script:

  SAMTOOLS
  MIN_MAPQ
  WINDOW_SIZE
  WINDOW_STEP
  MIN_REGION_MEAN_DEPTH
  DEPTH_RATIO_CUTOFF
  HET_DEL_RATIO_LOW
  HET_DEL_RATIO_HIGH
  MIN_CANDIDATE_WINDOWS
  MAX_MERGE_GAP_WINDOWS

Main output:

  depth_candidates.tsv

Additional output:

  *.window_depth.tsv
  *.region_mean_depth.tsv
  *.all_window_ratio.tsv

USAGE
}

