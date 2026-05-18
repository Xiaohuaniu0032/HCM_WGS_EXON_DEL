#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min max);

# -----------------------------
# Default parameters
# -----------------------------
my $ref_fa;
my $gene_txt;
my $exon_txt;
my $gene;
my $del_exon;

my $out_prefix  = "sim";
my $read_len    = 150;
my $depth       = 100;

my $insert_mean = 350;
my $insert_sd   = 50;

my $left_flank  = 200;
my $right_flank = 200;

my $error_rate  = 0.001;
my $seed        = 42;
my $samtools    = "samtools";

GetOptions(
    "ref=s"          => \$ref_fa,
    "gene-txt=s"     => \$gene_txt,
    "exon-txt=s"     => \$exon_txt,
    "gene=s"         => \$gene,
    "del-exon=s"     => \$del_exon,
    "out=s"          => \$out_prefix,
    "read-len=i"     => \$read_len,
    "depth=f"        => \$depth,
    "insert-mean=i"  => \$insert_mean,
    "insert-sd=i"    => \$insert_sd,
    "left-flank=i"   => \$left_flank,
    "right-flank=i"  => \$right_flank,
    "error-rate=f"   => \$error_rate,
    "seed=i"         => \$seed,
    "samtools=s"     => \$samtools,
) or die usage();

die usage() unless $ref_fa && $gene_txt && $exon_txt && $gene && $del_exon;

# -----------------------------
# Parameter checking
# -----------------------------
die "[ERROR] --read-len must be > 0\n" if $read_len <= 0;
die "[ERROR] --depth must be > 0\n" if $depth <= 0;
die "[ERROR] --insert-mean must be >= 2 * read_len\n"
    if $insert_mean < 2 * $read_len;
die "[ERROR] --insert-sd must be >= 0\n" if $insert_sd < 0;
die "[ERROR] --left-flank must be >= 0\n" if $left_flank < 0;
die "[ERROR] --right-flank must be >= 0\n" if $right_flank < 0;
die "[ERROR] --error-rate must be between 0 and 1\n"
    if $error_rate < 0 || $error_rate > 1;

srand($seed);

# -----------------------------
# 1. Parse gene annotation
# -----------------------------
my $gene_info = parse_gene_txt($gene_txt, $gene);

my ($transcript, $chrom, $gene_start, $gene_end, $strand, $exon_count) =
    @{$gene_info}{qw/transcript chrom start end strand exon_count/};

my $gene_region = "$chrom:$gene_start-$gene_end";

# -----------------------------
# 2. Parse deletion exon string
#    E15 or E12_14
# -----------------------------
my @target_exons = parse_del_exons($del_exon);

# -----------------------------
# 3. Parse all exons of target transcript
# -----------------------------
my @all_exon_records = parse_all_exons_of_transcript(
    $exon_txt,
    $gene,
    $transcript
);

die "[ERROR] No exon records found for $gene $transcript in $exon_txt\n"
    unless @all_exon_records;

# Select target exon records
my %target_exon_hash = map { $_ => 1 } @target_exons;

my @target_exon_records = grep {
    exists $target_exon_hash{ $_->{exon} }
} @all_exon_records;

die "[ERROR] No target exon records found for $gene $transcript $del_exon\n"
    unless @target_exon_records;

# Check if all requested exons exist
my %found_exons = map { $_->{exon} => 1 } @target_exon_records;

for my $e (@target_exons) {
    die "[ERROR] Requested exon E$e not found in $gene $transcript\n"
        unless exists $found_exons{$e};
}

# -----------------------------
# 4. Calculate deletion region with left/right flank
# -----------------------------
my (
    $del_start,
    $del_end,
    $raw_del_start,
    $raw_del_end,
    $left_limit,
    $right_limit
) = get_extended_deletion_region(
    all_exons    => \@all_exon_records,
    target_exons => \@target_exon_records,
    gene_start   => $gene_start,
    gene_end     => $gene_end,
    left_flank   => $left_flank,
    right_flank  => $right_flank,
);

my $del_len = $del_end - $del_start + 1;

die "[ERROR] Deletion region is outside gene region\n"
    if $del_start < $gene_start || $del_end > $gene_end;

# -----------------------------
# 5. Fetch gene sequence
# -----------------------------
my $wt_seq = fetch_fasta_region($ref_fa, $gene_region, $samtools);

$wt_seq = uc($wt_seq);
$wt_seq =~ s/[^ACGT]/N/g;

my $gene_len = length($wt_seq);

die "[ERROR] Fetched sequence length does not match gene interval length\n"
    if $gene_len != ($gene_end - $gene_start + 1);

# -----------------------------
# 6. Build deletion haplotype
# -----------------------------
my $del_offset_0 = $del_start - $gene_start;

my $del_seq =
    substr($wt_seq, 0, $del_offset_0) .
    substr($wt_seq, $del_offset_0 + $del_len);

my $del_hap_len = length($del_seq);

die "[ERROR] Read length is too long for target gene region\n"
    if $gene_len < $read_len * 2;

die "[ERROR] Read length is too long for deletion haplotype\n"
    if $del_hap_len < $read_len * 2;

# -----------------------------
# 7. Open output files
# -----------------------------
my $r1_fq = "$out_prefix.R1.fastq";
my $r2_fq = "$out_prefix.R2.fastq";
my $truth = "$out_prefix.truth.txt";

open my $R1, ">", $r1_fq or die "[ERROR] Cannot write $r1_fq: $!\n";
open my $R2, ">", $r2_fq or die "[ERROR] Cannot write $r2_fq: $!\n";
open my $TR, ">", $truth or die "[ERROR] Cannot write $truth: $!\n";

# -----------------------------
# 8. Write truth file
# -----------------------------
print $TR join("\t", qw/
Gene Transcript Chrom GeneStart GeneEnd Strand ExonCount
DelExon RawDelStart RawDelEnd ExtendedDelStart ExtendedDelEnd DelLen
LeftFlank RightFlank LeftBoundaryLimit RightBoundaryLimit
ReadLen Depth InsertMean InsertSD ErrorRate
/), "\n";

print $TR join("\t",
    $gene,
    $transcript,
    $chrom,
    $gene_start,
    $gene_end,
    $strand,
    $exon_count,
    $del_exon,
    $raw_del_start,
    $raw_del_end,
    $del_start,
    $del_end,
    $del_len,
    $left_flank,
    $right_flank,
    $left_limit,
    $right_limit,
    $read_len,
    $depth,
    $insert_mean,
    $insert_sd,
    $error_rate
), "\n";

# -----------------------------
# 9. Simulate heterozygous deletion
#    WT haplotype:  50%
#    DEL haplotype: 50%
# -----------------------------
my $wt_depth  = $depth / 2;
my $del_depth = $depth / 2;

my $wt_frag_n  = int(($wt_depth  * $gene_len)    / (2 * $read_len) + 0.5);
my $del_frag_n = int(($del_depth * $del_hap_len) / (2 * $read_len) + 0.5);

my $read_id = 1;

simulate_haplotype_reads(
    seq          => $wt_seq,
    hap_name     => "WT",
    frag_n       => $wt_frag_n,
    read_len     => $read_len,
    insert_mean  => $insert_mean,
    insert_sd    => $insert_sd,
    error_rate   => $error_rate,
    r1_fh        => $R1,
    r2_fh        => $R2,
    read_id_ref  => \$read_id,
    region_name  => $gene_region,
);

simulate_haplotype_reads(
    seq          => $del_seq,
    hap_name     => "DEL",
    frag_n       => $del_frag_n,
    read_len     => $read_len,
    insert_mean  => $insert_mean,
    insert_sd    => $insert_sd,
    error_rate   => $error_rate,
    r1_fh        => $R1,
    r2_fh        => $R2,
    read_id_ref  => \$read_id,
    region_name  => $gene_region,
);

close $R1;
close $R2;
close $TR;

# -----------------------------
# 10. Log information
# -----------------------------
print STDERR "[INFO] Gene: $gene\n";
print STDERR "[INFO] Transcript: $transcript\n";
print STDERR "[INFO] Strand: $strand\n";
print STDERR "[INFO] Gene region: $gene_region\n";
print STDERR "[INFO] Gene length: $gene_len bp\n";
print STDERR "[INFO] Deletion exon: $del_exon\n";
print STDERR "[INFO] Raw exon deletion region: $chrom:$raw_del_start-$raw_del_end\n";
print STDERR "[INFO] Left flank requested: $left_flank bp\n";
print STDERR "[INFO] Right flank requested: $right_flank bp\n";
print STDERR "[INFO] Left boundary limit: $left_limit\n";
print STDERR "[INFO] Right boundary limit: $right_limit\n";
print STDERR "[INFO] Extended deletion region: $chrom:$del_start-$del_end\n";
print STDERR "[INFO] Extended deletion length: $del_len bp\n";
print STDERR "[INFO] WT haplotype length: $gene_len bp\n";
print STDERR "[INFO] DEL haplotype length: $del_hap_len bp\n";
print STDERR "[INFO] Total depth: ${depth}X\n";
print STDERR "[INFO] WT depth: ${wt_depth}X\n";
print STDERR "[INFO] DEL depth: ${del_depth}X\n";
print STDERR "[INFO] WT fragments: $wt_frag_n\n";
print STDERR "[INFO] DEL fragments: $del_frag_n\n";
print STDERR "[INFO] Output R1: $r1_fq\n";
print STDERR "[INFO] Output R2: $r2_fq\n";
print STDERR "[INFO] Truth file: $truth\n";

exit 0;

# ============================================================
# Subroutines
# ============================================================

sub usage {
    return <<"USAGE";

Usage:
  perl simulate_gene_wgs_del.pl \\
    --ref hg19.fa \\
    --gene-txt RefSeq_MANE_Select.gene.txt \\
    --exon-txt RefSeq_MANE_Select.exon.txt \\
    --gene FHOD3 \\
    --del-exon E12_14 \\
    --left-flank 200 \\
    --right-flank 200 \\
    --out FHOD3_E12_14del \\
    --read-len 150 \\
    --depth 100

Required:
  --ref             Reference genome FASTA
  --gene-txt        RefSeq_MANE_Select.gene.txt
  --exon-txt        RefSeq_MANE_Select.exon.txt
  --gene            Target gene name, for example FHOD3
  --del-exon        Deleted exon, for example E15 or E12_14

Optional:
  --out             Output prefix, default: sim
  --read-len        PE read length, default: 150
  --depth           Total sequencing depth, default: 100
  --insert-mean     Insert size mean, default: 350
  --insert-sd       Insert size standard deviation, default: 50
  --left-flank      Left breakpoint flank size, default: 200
  --right-flank     Right breakpoint flank size, default: 200
  --error-rate      Sequencing error rate, default: 0.001
  --seed            Random seed, default: 42
  --samtools        samtools path, default: samtools

Output:
  prefix.R1.fastq
  prefix.R2.fastq
  prefix.truth.txt

Example:
  perl simulate_gene_wgs_del.pl \\
    --ref /ehpcdata/fulongfei/database/ref/hg19/hg19.fa \\
    --gene-txt RefSeq_MANE_Select.gene.txt \\
    --exon-txt RefSeq_MANE_Select.exon.txt \\
    --gene FHOD3 \\
    --del-exon E15 \\
    --left-flank 300 \\
    --right-flank 100 \\
    --out FHOD3_E15del \\
    --read-len 150 \\
    --depth 100

USAGE
}

sub parse_gene_txt {
    my ($file, $target_gene) = @_;

    open my $IN, "<", $file or die "[ERROR] Cannot open $file: $!\n";

    my $header = <$IN>;
    chomp $header;

    my @h = split /\t|\s+/, $header;

    my %idx;
    for my $i (0 .. $#h) {
        $idx{$h[$i]} = $i;
    }

    for my $col (qw/Gene Transcript Chrom Start End Strand ExonCount/) {
        die "[ERROR] Missing column $col in $file\n" unless exists $idx{$col};
    }

    while (<$IN>) {
        chomp;
        next if /^\s*$/;

        my @f = split /\t|\s+/;

        next unless $f[$idx{Gene}] eq $target_gene;

        close $IN;

        return {
            gene       => $f[$idx{Gene}],
            transcript => $f[$idx{Transcript}],
            chrom      => $f[$idx{Chrom}],
            start      => $f[$idx{Start}],
            end        => $f[$idx{End}],
            strand     => $f[$idx{Strand}],
            exon_count => $f[$idx{ExonCount}],
        };
    }

    close $IN;

    die "[ERROR] Gene $target_gene not found in $file\n";
}

sub parse_del_exons {
    my ($del) = @_;

    my $raw = $del;

    $del =~ s/^E//i;

    my @exons;

    if ($del =~ /^(\d+)$/) {
        @exons = ($1);
    }
    elsif ($del =~ /^(\d+)_(\d+)$/) {
        my ($s, $e) = ($1, $2);

        die "[ERROR] Invalid exon range: $raw\n" if $s > $e;

        @exons = ($s .. $e);
    }
    else {
        die "[ERROR] Invalid --del-exon format: $raw. Use E15 or E12_14\n";
    }

    return @exons;
}

sub parse_all_exons_of_transcript {
    my ($file, $target_gene, $target_tx) = @_;

    open my $IN, "<", $file or die "[ERROR] Cannot open $file: $!\n";

    my $header = <$IN>;
    chomp $header;

    my @h = split /\t|\s+/, $header;

    my %idx;
    for my $i (0 .. $#h) {
        $idx{$h[$i]} = $i;
    }

    for my $col (qw/Gene Transcript Exon Chrom Start End Strand/) {
        die "[ERROR] Missing column $col in $file\n" unless exists $idx{$col};
    }

    my @records;

    while (<$IN>) {
        chomp;
        next if /^\s*$/;

        my @f = split /\t|\s+/;

        next unless $f[$idx{Gene}] eq $target_gene;
        next unless $f[$idx{Transcript}] eq $target_tx;

        push @records, {
            gene       => $f[$idx{Gene}],
            transcript => $f[$idx{Transcript}],
            exon       => $f[$idx{Exon}],
            chrom      => $f[$idx{Chrom}],
            start      => $f[$idx{Start}],
            end        => $f[$idx{End}],
            strand     => $f[$idx{Strand}],
        };
    }

    close $IN;

    @records = sort {
        $a->{start} <=> $b->{start}
        ||
        $a->{end} <=> $b->{end}
    } @records;

    return @records;
}

sub get_extended_deletion_region {
    my %args = @_;

    my $all_exons_ref    = $args{all_exons};
    my $target_exons_ref = $args{target_exons};
    my $gene_start       = $args{gene_start};
    my $gene_end         = $args{gene_end};
    my $left_flank       = $args{left_flank};
    my $right_flank      = $args{right_flank};

    my $raw_start = min(map { $_->{start} } @$target_exons_ref);
    my $raw_end   = max(map { $_->{end}   } @$target_exons_ref);

    my $left_limit  = $gene_start;
    my $right_limit = $gene_end;

    for my $exon (@$all_exons_ref) {
        my $s = $exon->{start};
        my $e = $exon->{end};

        # Exons completely inside the raw deletion interval are target-region exons
        # and should not restrict the breakpoint.
        next if $s >= $raw_start && $e <= $raw_end;

        # The closest upstream exon restricts the left breakpoint.
        # The deletion start cannot be smaller than upstream_exon_end + 1.
        if ($e < $raw_start) {
            $left_limit = $e + 1 if $e + 1 > $left_limit;
        }

        # The closest downstream exon restricts the right breakpoint.
        # The deletion end cannot be larger than downstream_exon_start - 1.
        if ($s > $raw_end) {
            $right_limit = $s - 1 if $s - 1 < $right_limit;
        }
    }

    my $extended_start = $raw_start - $left_flank;
    my $extended_end   = $raw_end   + $right_flank;

    $extended_start = $left_limit  if $extended_start < $left_limit;
    $extended_end   = $right_limit if $extended_end   > $right_limit;

    die "[ERROR] Invalid extended deletion region: $extended_start-$extended_end\n"
        if $extended_start > $extended_end;

    return (
        $extended_start,
        $extended_end,
        $raw_start,
        $raw_end,
        $left_limit,
        $right_limit
    );
}

sub fetch_fasta_region {
    my ($fa, $region, $samtools) = @_;

    my $cmd = "$samtools faidx $fa $region";

    open my $FA, "-|", $cmd or die "[ERROR] Failed to run: $cmd\n";

    my $seq = "";

    while (<$FA>) {
        chomp;
        next if /^>/;
        $seq .= $_;
    }

    close $FA;

    die "[ERROR] No sequence fetched for $region. Please check FASTA index and chromosome name.\n"
        unless $seq;

    return $seq;
}

sub simulate_haplotype_reads {
    my %args = @_;

    my $seq         = $args{seq};
    my $hap_name    = $args{hap_name};
    my $frag_n      = $args{frag_n};
    my $read_len    = $args{read_len};
    my $insert_mean = $args{insert_mean};
    my $insert_sd   = $args{insert_sd};
    my $error_rate  = $args{error_rate};
    my $R1          = $args{r1_fh};
    my $R2          = $args{r2_fh};
    my $read_id_ref = $args{read_id_ref};
    my $region_name = $args{region_name};

    my $seq_len = length($seq);

    for my $i (1 .. $frag_n) {
        my $insert = sample_insert_size(
            $insert_mean,
            $insert_sd,
            $read_len,
            $seq_len
        );

        my $max_start = $seq_len - $insert + 1;

        next if $max_start < 1;

        my $frag_start_1 = int(rand($max_start)) + 1;
        my $frag_start_0 = $frag_start_1 - 1;

        my $fragment = substr($seq, $frag_start_0, $insert);

        my $r1_seq = substr($fragment, 0, $read_len);
        my $r2_seq = substr($fragment, $insert - $read_len, $read_len);

        $r2_seq = revcomp($r2_seq);

        $r1_seq = add_errors($r1_seq, $error_rate);
        $r2_seq = add_errors($r2_seq, $error_rate);

        my $qual = "I" x $read_len;

        my $id = $$read_id_ref;
        $$read_id_ref++;

        my $header =
            "\@SIM:${id}:hap=${hap_name}:region=${region_name}:start=${frag_start_1}:insert=${insert}";

        print $R1 "$header/1\n$r1_seq\n+\n$qual\n";
        print $R2 "$header/2\n$r2_seq\n+\n$qual\n";
    }
}

sub sample_insert_size {
    my ($mean, $sd, $read_len, $seq_len) = @_;

    my $insert;

    for (1 .. 1000) {
        $insert = int(rand_normal($mean, $sd) + 0.5);

        last if $insert >= 2 * $read_len && $insert <= $seq_len;
    }

    $insert = 2 * $read_len if $insert < 2 * $read_len;
    $insert = $seq_len      if $insert > $seq_len;

    return $insert;
}

sub rand_normal {
    my ($mean, $sd) = @_;

    return $mean if $sd == 0;

    my $u1 = rand();
    my $u2 = rand();

    $u1 = 1e-10 if $u1 == 0;

    my $z = sqrt(-2 * log($u1)) * cos(2 * 3.141592653589793 * $u2);

    return $mean + $sd * $z;
}

sub revcomp {
    my ($seq) = @_;

    $seq = reverse($seq);
    $seq =~ tr/ACGTNacgtn/TGCANtgcan/;

    return uc($seq);
}

sub add_errors {
    my ($seq, $err) = @_;

    return $seq if $err <= 0;

    my @base = split //, $seq;
    my @nt = qw/A C G T/;

    for my $i (0 .. $#base) {
        next if rand() > $err;
        next if $base[$i] eq "N";

        my @alt = grep { $_ ne $base[$i] } @nt;

        $base[$i] = $alt[int(rand(@alt))];
    }

    return join("", @base);
}

