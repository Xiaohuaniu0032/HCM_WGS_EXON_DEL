#!/usr/bin/env perl
use strict;
use warnings;

# Usage:
# perl calc_depth_ratio.pl depth.txt depth.ratio.tsv depth.summary.tsv

my ($infile, $outfile, $summary) = @ARGV;

die "Usage: perl $0 <depth.txt> <depth.ratio.tsv> <depth.summary.tsv>\n"
    unless @ARGV == 3;

open my $IN, "<", $infile or die "Cannot open $infile: $!\n";

my @depths;
my $sum = 0;
my $n   = 0;

while (my $line = <$IN>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    next if $line =~ /^#/;

    my @F = split /\s+/, $line;
    next unless @F >= 3;

    my $depth = $F[2];

    unless ($depth =~ /^-?\d+(\.\d+)?$/) {
        warn "Skip non-numeric depth line: $line\n";
        next;
    }

    push @depths, $depth;
    $sum += $depth;
    $n++;
}

close $IN;

die "No valid depth records found in $infile\n" if $n == 0;

# 计算均值
my $mean_depth = $sum / $n;

# 计算中位值
my @sorted = sort { $a <=> $b } @depths;
my $median_depth;

if ($n % 2 == 1) {
    $median_depth = $sorted[int($n / 2)];
} else {
    my $mid1 = $sorted[$n / 2 - 1];
    my $mid2 = $sorted[$n / 2];
    $median_depth = ($mid1 + $mid2) / 2;
}

# 输出 summary
open my $SUM, ">", $summary or die "Cannot write $summary: $!\n";
print $SUM "Metric\tValue\n";
print $SUM "Total_Sites\t$n\n";
print $SUM "Median_Depth\t", sprintf("%.6f", $median_depth), "\n";
print $SUM "Mean_Depth\t", sprintf("%.6f", $mean_depth), "\n";
close $SUM;

# 第二遍读取文件，计算每个位点 ratio
open $IN, "<", $infile or die "Cannot open $infile: $!\n";
open my $OUT, ">", $outfile or die "Cannot write $outfile: $!\n";

print $OUT join("\t",
    "#CHROM",
    "POS",
    "Depth",
    "Median_Depth",
    "Mean_Depth",
    "Ratio_By_Median",
    "Ratio_By_Mean"
), "\n";

while (my $line = <$IN>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    next if $line =~ /^#/;

    my @F = split /\s+/, $line;
    next unless @F >= 3;

    my ($chr, $pos, $depth) = @F[0, 1, 2];

    next unless $depth =~ /^-?\d+(\.\d+)?$/;

    my $ratio_median = $median_depth == 0 ? "NA" : sprintf("%.6f", $depth / $median_depth);
    my $ratio_mean   = $mean_depth   == 0 ? "NA" : sprintf("%.6f", $depth / $mean_depth);

    print $OUT join("\t",
        $chr,
        $pos,
        $depth,
        sprintf("%.6f", $median_depth),
        sprintf("%.6f", $mean_depth),
        $ratio_median,
        $ratio_mean
    ), "\n";
}

close $IN;
close $OUT;

print STDERR "Done.\n";
print STDERR "Total sites   : $n\n";
print STDERR "Median depth  : ", sprintf("%.6f", $median_depth), "\n";
print STDERR "Mean depth    : ", sprintf("%.6f", $mean_depth), "\n";
print STDERR "Output ratio  : $outfile\n";
print STDERR "Output summary: $summary\n";

