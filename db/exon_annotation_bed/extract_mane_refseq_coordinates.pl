#!/usr/bin/env perl
# -*- coding: utf-8 -*-

use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);
use Cwd qw(abs_path);

# ============================================================
# extract_mane_exon_bed.pl
#
# Function:
#   Extract exon-level and gene-level genomic coordinates for
#   RefSeq MANE Select transcripts from a GTF file.
#
# Output:
#   1. RefSeq_MANE_Select.exon.txt
#   2. RefSeq_MANE_Select.gene.txt
#
# Coordinate:
#   Input GTF : 1-based closed interval
#   Output TXT: 1-based closed interval
#
# Exon TXT columns:
#   Gene
#   Transcript
#   Exon
#   Chrom
#   Start
#   End
#   Strand
#
# Gene TXT columns:
#   Gene
#   Transcript
#   Chrom
#   Start
#   End
#   Strand
#   ExonCount
#
# ============================================================

my ($gtf, $mane, $outdir, $exon_out, $gene_out, $unmatched, $keep_version, $help);

GetOptions(
    "gtf=s"          => \$gtf,
    "mane=s"         => \$mane,
    "outdir=s"       => \$outdir,
    "exon-out=s"     => \$exon_out,
    "gene-out=s"     => \$gene_out,
    "unmatched=s"    => \$unmatched,
    "keep-version!"  => \$keep_version,
    "help"           => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless $gtf && $mane;

$gtf  = abs_path($gtf)  or die "[ERROR] GTF file not found: $gtf\n";
$mane = abs_path($mane) or die "[ERROR] MANE file not found: $mane\n";

$outdir ||= ".";
make_path($outdir) unless -d $outdir;
$outdir = abs_path($outdir);

$exon_out  ||= "$outdir/RefSeq_MANE_Select.exon.txt";
$gene_out  ||= "$outdir/RefSeq_MANE_Select.gene.txt";
$unmatched ||= "$outdir/RefSeq_MANE_Select.unmatched.txt";

# ------------------------------------------------------------
# Step 1. Read MANE Select transcript table
# ------------------------------------------------------------

my %mane_tx_to_gene;
my %mane_gene_to_tx;
my %mane_tx_raw;

open my $MANE, "<", $mane or die "[ERROR] Cannot open MANE file: $mane\n";

my $mane_line_no = 0;
my $loaded_mane  = 0;

while (my $line = <$MANE>) {
    chomp $line;
    $line =~ s/\r$//;
    $mane_line_no++;

    next if $line =~ /^\s*$/;
    next if $line =~ /^\s*#/;

    my @F = split /\t/, $line;

    if ($mane_line_no == 1) {
        my $header = join("\t", @F);
        if ($header =~ /gene|name|symbol/i && $header =~ /RefSeq|transcript|MANE/i) {
            next;
        }
    }

    # Expected format:
    #   Gene    RefSeq_transcript    MANE_status
    #
    # Compatible previous format:
    #   name    RefSeq_prot          MANE_status
    my $gene = $F[0] // "";
    my $tx   = $F[1] // "";
    my $stat = $F[2] // "";

    $gene =~ s/^\s+|\s+$//g;
    $tx   =~ s/^\s+|\s+$//g;
    $stat =~ s/^\s+|\s+$//g;

    next if $gene eq "" || $tx eq "";

    if ($stat ne "" && $stat !~ /MANE\s*Select/i) {
        next;
    }

    my $tx_norm = normalize_refseq_id($tx, $keep_version);

    $mane_tx_to_gene{$tx_norm} = $gene;
    $mane_gene_to_tx{$gene}    = $tx_norm;
    $mane_tx_raw{$tx_norm}     = $tx;

    $loaded_mane++;
}

close $MANE;

die "[ERROR] No MANE Select transcripts loaded from: $mane\n"
    if $loaded_mane == 0;

# ------------------------------------------------------------
# Step 2. Parse GTF exon records
# ------------------------------------------------------------

my @exon_records;
my %gene_range;
my %matched_genes;
my %matched_txs;
my %skipped_chr;

my $total_gtf_records          = 0;
my $total_gtf_exons            = 0;
my $noncanonical_chr_skipped   = 0;
my $matched_exons              = 0;
my $gene_name_mismatch_skipped = 0;

open my $GTF, "<", $gtf or die "[ERROR] Cannot open GTF file: $gtf\n";

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

    unless (is_canonical_chr($chr)) {
        $noncanonical_chr_skipped++;
        $skipped_chr{$chr} = 1;
        next;
    }

    my %attr = parse_gtf_attr($attr);

    my $gene_name     = $attr{"gene_name"}     // $attr{"gene"}      // "";
    my $transcript_id = $attr{"transcript_id"} // "";
    my $exon_number   = $attr{"exon_number"}   // ".";

    next if $transcript_id eq "";

    my $tx_norm = normalize_refseq_id($transcript_id, $keep_version);

    next unless exists $mane_tx_to_gene{$tx_norm};

    my $mane_gene = $mane_tx_to_gene{$tx_norm};

    if ($gene_name ne "" && $gene_name ne $mane_gene) {
        $gene_name_mismatch_skipped++;
        next;
    }

    my $out_gene = $gene_name ne "" ? $gene_name : $mane_gene;

    # GTF coordinates are already 1-based closed.
    my $out_start = $start;
    my $out_end   = $end;

    next if $out_start < 1;
    next if $out_end < $out_start;

    my $exon_label = format_exon_number($exon_number);

    push @exon_records, {
        gene       => $out_gene,
        transcript => $tx_norm,
        exon       => $exon_label,
        chr        => $chr,
        start      => $out_start,
        end        => $out_end,
        strand     => $strand,
    };

    $matched_exons++;
    $matched_genes{$out_gene} = 1;
    $matched_txs{$tx_norm}    = 1;

    my $key = join("\t", $out_gene, $tx_norm, $chr, $strand);

    if (!exists $gene_range{$key}) {
        $gene_range{$key} = {
            gene       => $out_gene,
            transcript => $tx_norm,
            chr        => $chr,
            start      => $out_start,
            end        => $out_end,
            strand     => $strand,
            exon_count => 1,
        };
    } else {
        $gene_range{$key}->{start} = $out_start if $out_start < $gene_range{$key}->{start};
        $gene_range{$key}->{end}   = $out_end   if $out_end   > $gene_range{$key}->{end};
        $gene_range{$key}->{exon_count}++;
    }
}

close $GTF;

# ------------------------------------------------------------
# Step 3. Write exon TXT
# ------------------------------------------------------------

open my $EXON, ">", $exon_out or die "[ERROR] Cannot write exon TXT: $exon_out\n";

print $EXON join("\t",
    qw(Gene Transcript Exon Chrom Start End Strand)
), "\n";

foreach my $r (sort exon_sort @exon_records) {
    print $EXON join("\t",
        $r->{gene},
        $r->{transcript},
        $r->{exon},
        $r->{chr},
        $r->{start},
        $r->{end},
        $r->{strand},
    ), "\n";
}

close $EXON;

# ------------------------------------------------------------
# Step 4. Write gene TXT
# ------------------------------------------------------------

open my $GENE, ">", $gene_out or die "[ERROR] Cannot write gene TXT: $gene_out\n";

print $GENE join("\t",
    qw(Gene Transcript Chrom Start End Strand ExonCount)
), "\n";

foreach my $key (
    sort {
        chr_rank($gene_range{$a}->{chr}) <=> chr_rank($gene_range{$b}->{chr})
        ||
        $gene_range{$a}->{start} <=> $gene_range{$b}->{start}
        ||
        $gene_range{$a}->{end} <=> $gene_range{$b}->{end}
        ||
        $gene_range{$a}->{gene} cmp $gene_range{$b}->{gene}
    } keys %gene_range
) {
    my $r = $gene_range{$key};

    print $GENE join("\t",
        $r->{gene},
        $r->{transcript},
        $r->{chr},
        $r->{start},
        $r->{end},
        $r->{strand},
        $r->{exon_count},
    ), "\n";
}

close $GENE;

# ------------------------------------------------------------
# Step 5. Write unmatched report
# ------------------------------------------------------------

open my $UN, ">", $unmatched or die "[ERROR] Cannot write unmatched file: $unmatched\n";

print $UN join("\t",
    qw(Gene RefSeqTranscript RawRefSeqTranscript Reason)
), "\n";

foreach my $gene (sort keys %mane_gene_to_tx) {
    my $tx = $mane_gene_to_tx{$gene};

    unless (exists $matched_txs{$tx}) {
        print $UN join("\t",
            $gene,
            $tx,
            $mane_tx_raw{$tx} // ".",
            "Not_found_in_GTF_exon_records"
        ), "\n";
    }
}

close $UN;

# ------------------------------------------------------------
# Step 6. Summary
# ------------------------------------------------------------

print STDERR "Done.\n";
print STDERR "MANE Select transcripts loaded: $loaded_mane\n";
print STDERR "Total GTF records: $total_gtf_records\n";
print STDERR "Total GTF exon records: $total_gtf_exons\n";
print STDERR "Non-canonical chromosome exon records skipped: $noncanonical_chr_skipped\n";
print STDERR "Gene-name mismatch exon records skipped: $gene_name_mismatch_skipped\n";
print STDERR "Matched MANE exon records: $matched_exons\n";
print STDERR "Matched genes: ", scalar(keys %matched_genes), "\n";
print STDERR "Matched transcripts: ", scalar(keys %matched_txs), "\n";
print STDERR "Output exon TXT: $exon_out\n";
print STDERR "Output gene TXT: $gene_out\n";
print STDERR "Unmatched MANE transcript report: $unmatched\n";

if ($matched_exons == 0) {
    print STDERR "\n[WARNING] No exon records were matched.\n";
    print STDERR "Possible reasons:\n";
    print STDERR "  1. RefSeq transcript IDs in GTF and MANE file are inconsistent.\n";
    print STDERR "  2. The GTF file does not contain NM_/NR_ transcript IDs.\n";
    print STDERR "  3. The MANE file is not tab-delimited.\n";
    print STDERR "  4. Transcript version handling is inconsistent.\n";
    print STDERR "     Try using --keep-version or removing --keep-version.\n";
    print STDERR "\n";
}

exit 0;

# ============================================================
# Subroutines
# ============================================================

sub is_canonical_chr {
    my ($chr) = @_;

    return 1 if $chr =~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/;
    return 1 if $chr =~ /^([1-9]|1[0-9]|2[0-2]|X|Y)$/;

    return 0;
}

sub normalize_refseq_id {
    my ($tx, $keep_version) = @_;

    $tx //= "";
    $tx =~ s/^\s+|\s+$//g;

    return $tx if $keep_version;

    # Examples:
    #   NM_001281740.3 -> NM_001281740
    #   NM_001509.3_2  -> NM_001509
    $tx =~ s/\.\d+_\d+$//;
    $tx =~ s/\.\d+$//;

    return $tx;
}

sub parse_gtf_attr {
    my ($attr_string) = @_;

    my %attr;

    # Standard GTF:
    # gene_id "xxx"; transcript_id "xxx"; exon_number "1";
    while ($attr_string =~ /(\S+)\s+"([^"]*)"/g) {
        $attr{$1} = $2;
    }

    # Also support key=value style.
    while ($attr_string =~ /(\S+)=([^;\s]+)/g) {
        $attr{$1} = $2 unless exists $attr{$1};
    }

    return %attr;
}

sub format_exon_number {
    my ($exon_number) = @_;

    $exon_number //= ".";
    $exon_number =~ s/^\s+|\s+$//g;

    return "." if $exon_number eq "" || $exon_number eq ".";

    return $exon_number;
}

sub chr_rank {
    my ($chr) = @_;

    $chr =~ s/^chr//i;

    return $chr if $chr =~ /^\d+$/;
    return 23 if uc($chr) eq "X";
    return 24 if uc($chr) eq "Y";
    return 25 if uc($chr) eq "M" || uc($chr) eq "MT";

    return 1000;
}

sub exon_sort {
    return chr_rank($a->{chr}) <=> chr_rank($b->{chr})
        || $a->{start} <=> $b->{start}
        || $a->{end} <=> $b->{end}
        || $a->{gene} cmp $b->{gene}
        || $a->{transcript} cmp $b->{transcript}
        || exon_rank($a->{exon}) <=> exon_rank($b->{exon});
}

sub exon_rank {
    my ($exon) = @_;
    return $exon if defined $exon && $exon =~ /^\d+$/;
    return 999999;
}

sub usage {
    return <<"USAGE";

Usage:
  perl extract_mane_exon_bed.pl \\
      --gtf hg19.ncbiRefSeq.gtf \\
      --mane RefSeq_MANE_Select.xls \\
      --outdir ./

Output:
  RefSeq_MANE_Select.exon.txt
  RefSeq_MANE_Select.gene.txt
  RefSeq_MANE_Select.unmatched.txt

Required arguments:
  --gtf          Input GTF file.
  --mane         MANE Select transcript table.

Optional arguments:
  --outdir       Output directory.
                 Default: current directory.

  --exon-out     Output exon TXT file.
                 Default: <outdir>/RefSeq_MANE_Select.exon.txt

  --gene-out     Output gene TXT file.
                 Default: <outdir>/RefSeq_MANE_Select.gene.txt

  --unmatched    Output unmatched transcript report.
                 Default: <outdir>/RefSeq_MANE_Select.unmatched.txt

  --keep-version Keep RefSeq transcript version.
                 Default: remove transcript version.
                 Example:
                   NM_001281740.3 -> NM_001281740

  --help         Show this help message.

Input MANE file format:
  Gene    RefSeq_transcript    MANE_status
  A1BG    NM_130786            MANE Select
  A2M     NM_000014            MANE Select

Exon TXT format:
  Gene    Transcript    Exon    Chrom    Start    End    Strand

Gene TXT format:
  Gene    Transcript    Chrom   Start    End      Strand   ExonCount

Coordinate:
  GTF input : 1-based closed interval
  TXT output: 1-based closed interval

Example:
  perl extract_mane_exon_bed.pl \\
      --gtf hg19.ncbiRefSeq.gtf \\
      --mane ../canonical_transcripts/RefSeq_MANE_Select.xls \\
      --outdir ./ \\
      2> extract_mane_exon_bed.log

USAGE
}

