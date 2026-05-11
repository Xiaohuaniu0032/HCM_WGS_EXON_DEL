#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my ($gtf, $mane, $out, $unmatched, $help);

GetOptions(
    "gtf=s"       => \$gtf,
    "mane=s"      => \$mane,
    "out=s"       => \$out,
    "unmatched=s" => \$unmatched,
    "help"        => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless ($gtf && $mane && $out);

# ------------------------------------------------------------
# Step 1. Read MANE Select transcript file
# ------------------------------------------------------------
my %mane_tx_to_gene;
my %mane_gene_to_tx;

open my $MANE, "<", $mane or die "[ERROR] Cannot open MANE file: $mane\n";

my $line_no = 0;

while (my $line = <$MANE>) {
    chomp $line;
    $line =~ s/\r$//;
    $line_no++;

    next if $line =~ /^\s*$/;

    my @F = split /\t/, $line;

    # Skip header
    if ($line_no == 1 && $F[0] =~ /^name$/i) {
        next;
    }

    my $gene = $F[0] // "";
    my $tx   = $F[1] // "";
    my $stat = $F[2] // "";

    $gene =~ s/^\s+|\s+$//g;
    $tx   =~ s/^\s+|\s+$//g;
    $stat =~ s/^\s+|\s+$//g;

    next if $gene eq "" || $tx eq "";

    # 只保留 MANE Select
    if ($stat ne "" && $stat !~ /MANE\s+Select/i) {
        next;
    }

    my $tx_norm = normalize_refseq_id($tx);

    $mane_tx_to_gene{$tx_norm} = $gene;
    $mane_gene_to_tx{$gene} = $tx_norm;
}

close $MANE;

die "[ERROR] No MANE transcripts loaded from $mane\n"
    if scalar(keys %mane_tx_to_gene) == 0;

# ------------------------------------------------------------
# Step 2. Parse GTF and extract exon records
# ------------------------------------------------------------
open my $GTF, "<", $gtf or die "[ERROR] Cannot open GTF file: $gtf\n";
open my $OUT, ">", $out or die "[ERROR] Cannot write output file: $out\n";

my $total_gtf_records        = 0;
my $total_gtf_exons          = 0;
my $noncanonical_chr_skipped = 0;
my $matched_exons            = 0;

my %matched_genes;
my %matched_txs;
my %skipped_chr;

while (my $line = <$GTF>) {
    chomp $line;
    $line =~ s/\r$//;

    next if $line =~ /^\s*$/;
    next if $line =~ /^#/;

    $total_gtf_records++;

    my @F = split /\t/, $line;
    next unless @F >= 9;

    my ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $attr) = @F;

    next unless $feature eq "exon";
    $total_gtf_exons++;

    # 只保留常规染色体 chr1-chr22, chrX, chrY
    unless (is_canonical_chr($chr)) {
        $noncanonical_chr_skipped++;
        $skipped_chr{$chr} = 1;
        next;
    }

    my %attr = parse_gtf_attr($attr);

    my $gene_name     = $attr{"gene_name"}     // "";
    my $transcript_id = $attr{"transcript_id"} // "";
    my $exon_number   = $attr{"exon_number"}   // ".";

    next if $transcript_id eq "";

    my $tx_norm = normalize_refseq_id($transcript_id);

    next unless exists $mane_tx_to_gene{$tx_norm};

    my $mane_gene = $mane_tx_to_gene{$tx_norm};

    # 建议加 gene_name 校验，避免同一个 RefSeq ID 异常映射到其他 gene
    next if $gene_name ne "" && $gene_name ne $mane_gene;

    my $out_gene = $gene_name ne "" ? $gene_name : $mane_gene;

    my $bed_start = $start - 1;
    my $bed_end   = $end;

    next if $bed_start < 0;
    next if $bed_end <= $bed_start;

    my $bed_name = join("|",
        $out_gene,
        $tx_norm,
        "exon" . $exon_number
    );

    print $OUT join("\t",
        $chr,
        $bed_start,
        $bed_end,
        $bed_name,
        ".",
        $strand
    ), "\n";

    $matched_exons++;
    $matched_genes{$out_gene} = 1;
    $matched_txs{$tx_norm} = 1;
}

close $GTF;
close $OUT;

# ------------------------------------------------------------
# Step 3. Optional unmatched MANE transcript report
# ------------------------------------------------------------
if ($unmatched) {
    open my $UN, ">", $unmatched or die "[ERROR] Cannot write unmatched file: $unmatched\n";

    print $UN join("\t", qw(Gene RefSeq_transcript Reason)), "\n";

    foreach my $gene (sort keys %mane_gene_to_tx) {
        my $tx = $mane_gene_to_tx{$gene};

        unless (exists $matched_txs{$tx}) {
            print $UN join("\t", $gene, $tx, "Not_found_in_canonical_GTF_exon_records"), "\n";
        }
    }

    close $UN;
}

# ------------------------------------------------------------
# Step 4. Summary
# ------------------------------------------------------------
print STDERR "Done.\n";
print STDERR "MANE transcripts loaded: ", scalar(keys %mane_tx_to_gene), "\n";
print STDERR "Total GTF records: $total_gtf_records\n";
print STDERR "Total GTF exon records: $total_gtf_exons\n";
print STDERR "Non-canonical chromosome exon records skipped: $noncanonical_chr_skipped\n";
print STDERR "Matched MANE exon records: $matched_exons\n";
print STDERR "Matched genes: ", scalar(keys %matched_genes), "\n";
print STDERR "Matched transcripts: ", scalar(keys %matched_txs), "\n";
print STDERR "Output BED: $out\n";

if ($unmatched) {
    print STDERR "Unmatched MANE transcript report: $unmatched\n";
}

if ($matched_exons == 0) {
    print STDERR "\n[WARNING] No exon records were matched.\n";
    print STDERR "Possible reasons:\n";
    print STDERR "  1. RefSeq transcript IDs in GTF and MANE file are inconsistent.\n";
    print STDERR "  2. The GTF file does not contain NM_/NR_ transcript IDs.\n";
    print STDERR "  3. The MANE file is not tab-delimited.\n";
    print STDERR "  4. All matched transcripts are located on non-canonical chromosomes.\n";
    print STDERR "\n";
}

# ============================================================
# Subroutines
# ============================================================

sub is_canonical_chr {
    my ($chr) = @_;

    return 1 if $chr =~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/;

    return 0;
}

sub normalize_refseq_id {
    my ($tx) = @_;

    $tx =~ s/^\s+|\s+$//g;

    # 处理类似 NM_001509.3_2 的情况
    # NM_001509.3_2 -> NM_001509
    $tx =~ s/\.\d+_\d+$//;

    # 处理普通 RefSeq 版本号
    # NM_130786.4 -> NM_130786
    $tx =~ s/\.\d+$//;

    return $tx;
}

sub parse_gtf_attr {
    my ($attr_string) = @_;

    my %attr;

    while ($attr_string =~ /(\S+)\s+"([^"]*)"/g) {
        $attr{$1} = $2;
    }

    return %attr;
}

sub usage {
    return <<"USAGE";

Usage:
    perl extract_mane_exon_bed.pl \\
        --gtf hg19.ncbiRefSeq.gtf \\
        --mane RefSeq_MANE_Select.xls \\
        --out RefSeq_MANE_Select.exon.bed

Optional:
    --unmatched RefSeq_MANE_Select.unmatched.tsv

Input MANE file format:
    name    RefSeq_prot    MANE_status
    A1BG    NM_130786      MANE Select
    A2M     NM_000014      MANE Select

Output BED6 format:
    chrom    start0    end    gene|transcript|exon_number    .    strand

Only canonical chromosomes are retained:
    chr1-chr22, chrX, chrY

Coordinate:
    GTF: 1-based closed
    BED: 0-based half-open

Example:
    perl extract_mane_exon_bed.pl \\
        --gtf hg19.ncbiRefSeq.gtf \\
        --mane ../canonical_transcripts/RefSeq_MANE_Select.xls \\
        --out RefSeq_MANE_Select.exon.bed \\
        --unmatched RefSeq_MANE_Select.unmatched.tsv \\
        2> log

USAGE
}

