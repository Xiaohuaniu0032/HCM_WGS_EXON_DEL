#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use File::Spec;

# ============================================================
# gene_exon_intron_info.pl
#
# Purpose:
#   Given a gene name, extract all exons of this gene from
#   RefSeq_MANE_Select.exon.txt and calculate exon length,
#   genomic left/right intron length, and transcript-direction
#   upstream/downstream intron length.
#
# Definitions:
#
#   Exon_Length = Exon_End - Exon_Start + 1
#
#   Genomic_Left_Intron_Length:
#       Current exon start - previous genomic exon end - 1
#
#   Genomic_Right_Intron_Length:
#       Next genomic exon start - current exon end - 1
#
#   For + strand:
#       Transcript_Upstream_Intron_Length   = Genomic_Left_Intron_Length
#       Transcript_Downstream_Intron_Length = Genomic_Right_Intron_Length
#
#   For - strand:
#       Transcript_Upstream_Intron_Length   = Genomic_Right_Intron_Length
#       Transcript_Downstream_Intron_Length = Genomic_Left_Intron_Length
#
# Notes:
#   - Genomic left/right are always defined by chromosome coordinate order.
#   - Upstream/downstream are defined by transcript direction and therefore
#     depend on strand.
#   - The first/last exon will have NA for one side of the intron.
#
# Usage:
#   perl gene_exon_intron_info.pl \
#       --gene FHOD3 \
#       --exon-txt db/exon_annotation_bed/RefSeq_MANE_Select.exon.txt \
#       --outfile FHOD3.exon_intron_info.tsv
#
# Minimal usage if this script is placed in bin/:
#   perl bin/gene_exon_intron_info.pl --gene FHOD3
#
# ============================================================

my $gene;
my $exon_txt;
my $outfile;
my $help;

GetOptions(
    "gene=s"     => \$gene,
    "exon-txt=s" => \$exon_txt,
    "outfile=s"  => \$outfile,
    "help"       => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless defined $gene && $gene ne "";

# ------------------------------------------------------------
# Default annotation file
# Assume this script is placed in bin/
# ------------------------------------------------------------
if (!defined $exon_txt || $exon_txt eq "") {
    $exon_txt = File::Spec->catfile(
        $Bin,
        "..",
        "db",
        "exon_annotation_bed",
        "RefSeq_MANE_Select.exon.txt"
    );
}

die "[ERROR] Exon annotation file not found or empty: $exon_txt\n"
    unless -s $exon_txt;

# ------------------------------------------------------------
# Output handle
# ------------------------------------------------------------
my $OUT;
if (defined $outfile && $outfile ne "") {
    open $OUT, ">", $outfile
        or die "[ERROR] Cannot write outfile: $outfile\n";
}
else {
    $OUT = *STDOUT;
}

# ------------------------------------------------------------
# Read exon annotation file
# ------------------------------------------------------------
open my $IN, "<", $exon_txt
    or die "[ERROR] Cannot open exon annotation file: $exon_txt\n";

my $header = <$IN>;
die "[ERROR] Empty exon annotation file: $exon_txt\n" unless defined $header;

chomp $header;
$header =~ s/\r$//;

my @headers = split /\t/, $header;
my %col;

for my $i (0 .. $#headers) {
    my $name = $headers[$i];
    $name =~ s/^\s+|\s+$//g;
    $col{$name} = $i;
}

# ------------------------------------------------------------
# Flexible column detection
# Compatible with common column names
# ------------------------------------------------------------
my $gene_col = find_col(
    \%col,
    qw(
        Gene
        Gene_Name
        GeneSymbol
        Gene_Symbol
        gene
        gene_name
        gene_symbol
    )
);

my $transcript_col = find_col(
    \%col,
    qw(
        Transcript
        Transcript_ID
        TranscriptId
        Transcript_Name
        transcript
        transcript_id
        transcript_name
    )
);

my $chrom_col = find_col(
    \%col,
    qw(
        Chrom
        Chr
        chromosome
        chrom
        chr
    )
);

my $start_col = find_col(
    \%col,
    qw(
        Start
        Exon_Start
        exon_start
        start
    )
);

my $end_col = find_col(
    \%col,
    qw(
        End
        Exon_End
        exon_end
        end
    )
);

my $strand_col = find_col_optional(
    \%col,
    qw(
        Strand
        strand
    )
);

my $exon_col = find_col_optional(
    \%col,
    qw(
        Exon
        Exon_ID
        Exon_Name
        exon
        exon_id
        exon_name
    )
);

my $exon_no_col = find_col_optional(
    \%col,
    qw(
        Exon_Number
        Exon_No
        Exon_Order
        exon_number
        exon_no
        exon_order
    )
);

die "[ERROR] Required column for gene not found in header:\n$header\n"
    unless defined $gene_col;

die "[ERROR] Required column for transcript not found in header:\n$header\n"
    unless defined $transcript_col;

die "[ERROR] Required column for chrom not found in header:\n$header\n"
    unless defined $chrom_col;

die "[ERROR] Required column for start not found in header:\n$header\n"
    unless defined $start_col;

die "[ERROR] Required column for end not found in header:\n$header\n"
    unless defined $end_col;

# ------------------------------------------------------------
# Store exon records by transcript
# ------------------------------------------------------------
my %by_transcript;

while (my $line = <$IN>) {
    chomp $line;
    $line =~ s/\r$//;

    next if $line =~ /^\s*$/;
    next if $line =~ /^#/;

    my @f = split /\t/, $line, -1;

    my $this_gene = get_field(\@f, $gene_col);
    next unless defined $this_gene;
    next unless uc($this_gene) eq uc($gene);

    my $transcript = get_field(\@f, $transcript_col);
    my $chrom      = get_field(\@f, $chrom_col);
    my $start      = get_field(\@f, $start_col);
    my $end        = get_field(\@f, $end_col);

    next unless defined $transcript && $transcript ne "";
    next unless defined $chrom      && $chrom      ne "";
    next unless defined $start     && $start     =~ /^\d+$/;
    next unless defined $end       && $end       =~ /^\d+$/;

    my $strand = "NA";
    if (defined $strand_col) {
        my $tmp = get_field(\@f, $strand_col);
        $strand = $tmp if defined $tmp && $tmp ne "";
    }

    my $exon_id = "NA";
    if (defined $exon_col) {
        my $tmp = get_field(\@f, $exon_col);
        $exon_id = $tmp if defined $tmp && $tmp ne "";
    }

    my $exon_number = "NA";
    if (defined $exon_no_col) {
        my $tmp = get_field(\@f, $exon_no_col);
        $exon_number = $tmp if defined $tmp && $tmp ne "";
    }

    # If Exon_ID is missing but Exon_Number exists
    if ($exon_id eq "NA" && $exon_number ne "NA") {
        $exon_id = "EX" . $exon_number;
    }

    # If Exon_Number is missing but Exon_ID contains EX12-like pattern
    if ($exon_number eq "NA" && $exon_id =~ /EX(\d+)/i) {
        $exon_number = $1;
    }

    # Normalize coordinate order just in case
    my $s = int($start);
    my $e = int($end);

    if ($s > $e) {
        my $tmp = $s;
        $s = $e;
        $e = $tmp;
    }

    push @{ $by_transcript{$transcript} }, {
        Gene        => $this_gene,
        Transcript  => $transcript,
        Chrom       => $chrom,
        Strand      => $strand,
        Exon_ID     => $exon_id,
        Exon_Number => $exon_number,
        Start       => $s,
        End         => $e,
    };
}

close $IN;

# ------------------------------------------------------------
# Output header
# ------------------------------------------------------------
print $OUT join("\t",
    "Gene",
    "Transcript",
    "Chrom",
    "Strand",
    "Exon_ID",
    "Exon_Number",
    "Exon_Start",
    "Exon_End",
    "Exon_Length",
    "Genomic_Left_Intron_Length",
    "Genomic_Right_Intron_Length",
    "Transcript_Upstream_Intron_Length",
    "Transcript_Downstream_Intron_Length",
    "Previous_Genomic_Exon_ID",
    "Previous_Genomic_Exon_End",
    "Next_Genomic_Exon_ID",
    "Next_Genomic_Exon_Start",
    "Upstream_Exon_ID",
    "Downstream_Exon_ID"
), "\n";

# ------------------------------------------------------------
# Calculate exon and intron information
# ------------------------------------------------------------
my $found = 0;

for my $transcript (sort keys %by_transcript) {

    my @exons = sort {
        $a->{Start} <=> $b->{Start}
            ||
        $a->{End} <=> $b->{End}
            ||
        $a->{Exon_ID} cmp $b->{Exon_ID}
    } @{ $by_transcript{$transcript} };

    next unless @exons;
    $found = 1;

    for my $i (0 .. $#exons) {

        my $e = $exons[$i];

        my $exon_len = $e->{End} - $e->{Start} + 1;

        my $genomic_left_intron_len  = "NA";
        my $genomic_right_intron_len = "NA";

        my $prev_genomic_exon_id  = "NA";
        my $prev_genomic_exon_end = "NA";

        my $next_genomic_exon_id    = "NA";
        my $next_genomic_exon_start = "NA";

        my $upstream_intron_len   = "NA";
        my $downstream_intron_len = "NA";

        my $upstream_exon_id   = "NA";
        my $downstream_exon_id = "NA";

        # --------------------------------------------------------
        # Genomic left side:
        # previous exon in increasing genomic coordinate order
        # --------------------------------------------------------
        if ($i > 0) {
            my $prev = $exons[$i - 1];

            $genomic_left_intron_len = $e->{Start} - $prev->{End} - 1;
            $genomic_left_intron_len = 0 if $genomic_left_intron_len < 0;

            $prev_genomic_exon_id  = $prev->{Exon_ID};
            $prev_genomic_exon_end = $prev->{End};
        }

        # --------------------------------------------------------
        # Genomic right side:
        # next exon in increasing genomic coordinate order
        # --------------------------------------------------------
        if ($i < $#exons) {
            my $next = $exons[$i + 1];

            $genomic_right_intron_len = $next->{Start} - $e->{End} - 1;
            $genomic_right_intron_len = 0 if $genomic_right_intron_len < 0;

            $next_genomic_exon_id    = $next->{Exon_ID};
            $next_genomic_exon_start = $next->{Start};
        }

        # --------------------------------------------------------
        # Transcript direction-aware upstream/downstream
        # --------------------------------------------------------
        if ($e->{Strand} eq "+") {

            $upstream_intron_len   = $genomic_left_intron_len;
            $downstream_intron_len = $genomic_right_intron_len;

            $upstream_exon_id   = $prev_genomic_exon_id;
            $downstream_exon_id = $next_genomic_exon_id;
        }
        elsif ($e->{Strand} eq "-") {

            $upstream_intron_len   = $genomic_right_intron_len;
            $downstream_intron_len = $genomic_left_intron_len;

            $upstream_exon_id   = $next_genomic_exon_id;
            $downstream_exon_id = $prev_genomic_exon_id;
        }
        else {

            # Unknown strand
            $upstream_intron_len   = "NA";
            $downstream_intron_len = "NA";

            $upstream_exon_id   = "NA";
            $downstream_exon_id = "NA";
        }

        print $OUT join("\t",
            $e->{Gene},
            $e->{Transcript},
            $e->{Chrom},
            $e->{Strand},
            $e->{Exon_ID},
            $e->{Exon_Number},
            $e->{Start},
            $e->{End},
            $exon_len,
            $genomic_left_intron_len,
            $genomic_right_intron_len,
            $upstream_intron_len,
            $downstream_intron_len,
            $prev_genomic_exon_id,
            $prev_genomic_exon_end,
            $next_genomic_exon_id,
            $next_genomic_exon_start,
            $upstream_exon_id,
            $downstream_exon_id
        ), "\n";
    }
}

close $OUT if defined $outfile && $outfile ne "";

if (!$found) {
    die "[ERROR] No exon records found for gene: $gene in $exon_txt\n";
}

# ============================================================
# Subroutines
# ============================================================

sub get_field {
    my ($array_ref, $idx) = @_;
    return undef unless defined $idx;
    return undef if $idx > $#$array_ref;

    my $v = $array_ref->[$idx];
    return undef unless defined $v;

    $v =~ s/^\s+|\s+$//g;
    return $v;
}

sub find_col {
    my ($col_ref, @names) = @_;

    for my $name (@names) {
        return $col_ref->{$name} if exists $col_ref->{$name};
    }

    return undef;
}

sub find_col_optional {
    my ($col_ref, @names) = @_;

    for my $name (@names) {
        return $col_ref->{$name} if exists $col_ref->{$name};
    }

    return undef;
}

sub usage {
    return <<"USAGE";

Usage:
  perl gene_exon_intron_info.pl --gene GENE_NAME [options]

Required:
  --gene        Gene symbol, for example FHOD3, MYBPC3, MYH7

Optional:
  --exon-txt    Exon annotation TXT file
                Default:
                ../db/exon_annotation_bed/RefSeq_MANE_Select.exon.txt

  --outfile     Output TSV file
                Default: STDOUT

  --help        Show help message

Example:
  perl bin/gene_exon_intron_info.pl \\
      --gene FHOD3 \\
      --exon-txt db/exon_annotation_bed/RefSeq_MANE_Select.exon.txt \\
      --outfile FHOD3.exon_intron_info.tsv

Minimal example:
  perl bin/gene_exon_intron_info.pl --gene FHOD3 > FHOD3.exon_intron_info.tsv

Output columns:
  Gene
  Transcript
  Chrom
  Strand
  Exon_ID
  Exon_Number
  Exon_Start
  Exon_End
  Exon_Length
  Genomic_Left_Intron_Length
  Genomic_Right_Intron_Length
  Transcript_Upstream_Intron_Length
  Transcript_Downstream_Intron_Length
  Previous_Genomic_Exon_ID
  Previous_Genomic_Exon_End
  Next_Genomic_Exon_ID
  Next_Genomic_Exon_Start
  Upstream_Exon_ID
  Downstream_Exon_ID

USAGE
}


