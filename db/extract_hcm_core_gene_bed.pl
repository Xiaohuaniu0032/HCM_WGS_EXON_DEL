#!/usr/bin/env perl
use strict;
use warnings;

# 用法检查
if (@ARGV != 3) {
    die "Usage: perl $0 hcm_core_genes.txt RefSeq_MANE_Select.exon.bed output.bed\n";
}

my ($gene_file, $bed_file, $out_file) = @ARGV;

# 读取 HCM core gene 列表
my %core_genes;

open my $GF, "<", $gene_file or die "Cannot open $gene_file: $!\n";

while (my $line = <$GF>) {
    chomp $line;
    next if $line =~ /^\s*$/;

    my @fields = split /\s+/, $line;

    # 跳过表头
    next if $fields[0] eq "Gene";

    my $gene = $fields[0];

    # 去掉可能存在的空格
    $gene =~ s/^\s+|\s+$//g;

    $core_genes{$gene} = 1;
}

close $GF;

# 从 BED 文件中提取对应基因
open my $BED, "<", $bed_file or die "Cannot open $bed_file: $!\n";
open my $OUT, ">", $out_file or die "Cannot write $out_file: $!\n";

while (my $line = <$BED>) {
    chomp $line;
    next if $line =~ /^\s*$/;

    my @fields = split /\t/, $line;

    # BED 至少需要4列
    next unless @fields >= 4;

    my $info = $fields[3];

    # 第4列格式示例：IL9R|NM_002186|exon1
    my ($gene, $transcript, $exon) = split /\|/, $info;

    next unless defined $gene;

    if (exists $core_genes{$gene}) {
        print $OUT $line, "\n";
    }
}

close $BED;
close $OUT;

print "Done. Output written to: $out_file\n";

