#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);
use POSIX qw(strftime);

# ============================================================
# Script: simulate_random_gene_wgs_del.pl
#
# Purpose:
#   Randomly generate intragenic deletions for each gene using
#   predefined deletion size bins, construct WT/DEL haplotypes,
#   and simulate heterozygous paired-end FASTQ reads.
#
# Input gene file format:
#   Gene    Transcript      Chrom   Start   End     Strand  ExonCount
#   OR4F5   NM_001005484    chr1    65419   71585   +       3
#
# Main outputs:
#   outdir/random_gene_del.all.truth.tsv
#   outdir/<Sample>/<Sample>.R1.fastq
#   outdir/<Sample>/<Sample>.R2.fastq
#   outdir/<Sample>/<Sample>.truth.txt
#
# Notes:
#   - WT haplotype contributes depth / 2
#   - DEL haplotype contributes depth / 2
#   - The simulated sample represents a heterozygous deletion
#   - Deletions are generated within gene coordinates only
#   - Exon-overlap status is not considered
# ============================================================

my $ref_fa;
my $gene_txt;
my $gene_list;
my $outdir = "random_gene_del_sim";
my $samtools = "samtools";

my $replicates_per_bin = 3;
my $seed = 20260527;
my $sample_prefix = "HCM";

my $read_len = 150;
my $depth = 30;
my $insert_mean = 350;
my $insert_sd = 50;
my $error_rate = 0.001;

my $max_genes = 0;

GetOptions(
    "ref=s"              => \$ref_fa,
    "gene-txt=s"         => \$gene_txt,
    "gene-list=s"        => \$gene_list,
    "outdir=s"           => \$outdir,
    "samtools=s"         => \$samtools,
    "replicates=i"       => \$replicates_per_bin,
    "seed=i"             => \$seed,
    "sample-prefix=s"    => \$sample_prefix,
    "read-len=i"         => \$read_len,
    "depth=f"            => \$depth,
    "insert-mean=f"      => \$insert_mean,
    "insert-sd=f"        => \$insert_sd,
    "error-rate=f"       => \$error_rate,
    "max-genes=i"        => \$max_genes,
) or die usage();

die usage() unless defined $ref_fa && defined $gene_txt;

die "[ERROR] Reference FASTA not found: $ref_fa\n" unless -e $ref_fa;
die "[ERROR] Gene TXT not found: $gene_txt\n" unless -e $gene_txt;

make_path($outdir) unless -d $outdir;

my $time = strftime("%Y-%m-%d %H:%M:%S", localtime);

# ------------------------------------------------------------
# Deletion size bins
# ------------------------------------------------------------
my @bins = (
    { name => "50_100",      min => 50,   max => 100 },
    { name => "100_200",     min => 100,  max => 200 },
    { name => "200_500",     min => 200,  max => 500 },
    { name => "500_1000",    min => 500,  max => 1000 },
    { name => "1000_2000",   min => 1000, max => 2000 },
    { name => "2000_5000",   min => 2000, max => 5000 },
    { name => "5000_10000",  min => 5000, max => 10000 },
);

# ------------------------------------------------------------
# Optional target gene list
# ------------------------------------------------------------
my %target_gene;
if (defined $gene_list) {
    open my $GL, "<", $gene_list
        or die "[ERROR] Cannot open gene list: $gene_list\n";

    while (my $line = <$GL>) {
        chomp $line;
        next if $line =~ /^\s*$/;
        next if $line =~ /^#/;

        my @f = split /\t/, $line;
        my $g = $f[0];
        $g =~ s/^\s+|\s+$//g;
        next if $g eq "" || $g =~ /^Gene$/i;

        $target_gene{$g} = 1;
    }
    close $GL;

    my $n = scalar keys %target_gene;
    die "[ERROR] No valid genes found in gene list: $gene_list\n" if $n == 0;
    print STDERR "[INFO] Target gene list loaded: $n genes\n";
}

# ------------------------------------------------------------
# Output all-truth table
# ------------------------------------------------------------
my $all_truth = "$outdir/random_gene_del.all.truth.tsv";
open my $ALL, ">", $all_truth
    or die "[ERROR] Cannot write all truth file: $all_truth\n";

print $ALL join("\t",
    "Sample",
    "Gene",
    "Transcript",
    "Chrom",
    "Strand",
    "Gene_Start",
    "Gene_End",
    "Gene_Length",
    "ExonCount",
    "Length_Bin",
    "Bin_Min_Size",
    "Bin_Max_Size",
    "Replicate",
    "Del_Start",
    "Del_End",
    "Del_Size_bp",
    "WT_Haplotype_Length",
    "DEL_Haplotype_Length",
    "Depth",
    "WT_Depth",
    "DEL_Depth",
    "Read_Length",
    "Insert_Mean",
    "Insert_SD",
    "Error_Rate",
    "WT_Fragment_Number",
    "DEL_Fragment_Number",
    "Seed",
    "R1_FASTQ",
    "R2_FASTQ"
), "\n";

# ------------------------------------------------------------
# Read gene coordinate file
# ------------------------------------------------------------
open my $IN, "<", $gene_txt
    or die "[ERROR] Cannot open gene file: $gene_txt\n";

my $header = <$IN>;
chomp $header;

my @header_cols = split /\t/, $header;
my %col;
for my $i (0 .. $#header_cols) {
    $header_cols[$i] =~ s/^\s+|\s+$//g;
    $col{$header_cols[$i]} = $i;
}

for my $required (qw/Gene Transcript Chrom Start End Strand ExonCount/) {
    die "[ERROR] Required column not found in gene file: $required\n"
        unless exists $col{$required};
}

my $gene_count = 0;
my $event_count = 0;
my $skip_count = 0;

print STDERR "[INFO] Simulation started: $time\n";
print STDERR "[INFO] Reference FASTA       : $ref_fa\n";
print STDERR "[INFO] Gene TXT             : $gene_txt\n";
print STDERR "[INFO] Output directory     : $outdir\n";
print STDERR "[INFO] Read length          : $read_len\n";
print STDERR "[INFO] Total depth          : $depth\n";
print STDERR "[INFO] Insert mean          : $insert_mean\n";
print STDERR "[INFO] Insert SD            : $insert_sd\n";
print STDERR "[INFO] Error rate           : $error_rate\n";
print STDERR "[INFO] Replicates per bin   : $replicates_per_bin\n";
print STDERR "[INFO] Random seed base     : $seed\n";

while (my $line = <$IN>) {
    chomp $line;
    next if $line =~ /^\s*$/;

    my @f = split /\t/, $line;

    my $gene       = clean_field($f[$col{"Gene"}]);
    my $transcript = clean_field($f[$col{"Transcript"}]);
    my $chrom      = clean_field($f[$col{"Chrom"}]);
    my $start      = clean_field($f[$col{"Start"}]);
    my $end        = clean_field($f[$col{"End"}]);
    my $strand     = clean_field($f[$col{"Strand"}]);
    my $exon_count = clean_field($f[$col{"ExonCount"}]);

    next if $gene eq "";
    next if defined $gene_list && !exists $target_gene{$gene};

    die "[ERROR] Invalid coordinate for gene $gene: Start=$start End=$end\n"
        unless $start =~ /^\d+$/ && $end =~ /^\d+$/;

    if ($start > $end) {
        my $tmp = $start;
        $start = $end;
        $end = $tmp;
    }

    my $gene_len = $end - $start + 1;

    if ($gene_len < $read_len) {
        warn "[WARN] Skip gene $gene: gene length $gene_len < read length $read_len\n";
        next;
    }

    $gene_count++;
    last if $max_genes > 0 && $gene_count > $max_genes;

    print STDERR "[INFO] Processing gene: $gene ($chrom:$start-$end, length=$gene_len)\n";

    my $wt_seq = fetch_fasta_sequence($samtools, $ref_fa, $chrom, $start, $end);
    my $wt_len = length($wt_seq);

    die "[ERROR] FASTA length mismatch for $gene: expected $gene_len, got $wt_len\n"
        if $wt_len != $gene_len;

    for my $bin (@bins) {
        my $bin_name = $bin->{name};
        my $min_size = $bin->{min};
        my $max_size = $bin->{max};

        if ($gene_len < $min_size) {
            warn "[WARN] Skip gene $gene bin $bin_name: gene length $gene_len < bin min size $min_size\n";
            $skip_count += $replicates_per_bin;
            next;
        }

        my $effective_max = $max_size;
        $effective_max = $gene_len if $gene_len < $max_size;

        for my $rep (1 .. $replicates_per_bin) {

            $event_count++;
            my $event_seed = $seed + $event_count;
            srand($event_seed);

            my $del_size = random_int($min_size, $effective_max);
            my $max_del_start = $end - $del_size + 1;

            my $del_start = random_int($start, $max_del_start);
            my $del_end   = $del_start + $del_size - 1;

            my $del_offset_0 = $del_start - $start;

            my $del_seq = substr($wt_seq, 0, $del_offset_0)
                        . substr($wt_seq, $del_offset_0 + $del_size);

            my $del_hap_len = length($del_seq);

            if ($del_hap_len < $read_len) {
                warn "[WARN] Skip event: $gene $bin_name rep$rep because DEL haplotype length $del_hap_len < read length $read_len\n";
                $skip_count++;
                next;
            }

            my $safe_gene = sanitize_name($gene);
            my $sample = join("_",
                $sample_prefix,
                $safe_gene,
                "DEL",
                $bin_name,
                "rep$rep"
            );

            my $sample_dir = "$outdir/$sample";
            make_path($sample_dir) unless -d $sample_dir;

            my $r1_fastq = "$sample_dir/$sample.R1.fastq";
            my $r2_fastq = "$sample_dir/$sample.R2.fastq";
            my $truth_txt = "$sample_dir/$sample.truth.txt";

            my $wt_depth  = $depth / 2;
            my $del_depth = $depth / 2;

            my $wt_frag_n  = int(($wt_depth  * $wt_len)      / (2 * $read_len) + 0.5);
            my $del_frag_n = int(($del_depth * $del_hap_len) / (2 * $read_len) + 0.5);

            open my $R1, ">", $r1_fastq
                or die "[ERROR] Cannot write R1 FASTQ: $r1_fastq\n";
            open my $R2, ">", $r2_fastq
                or die "[ERROR] Cannot write R2 FASTQ: $r2_fastq\n";

            simulate_haplotype_reads(
                seq          => $wt_seq,
                hap_name     => "WT",
                chrom        => $chrom,
                gene_start   => $start,
                gene_end     => $end,
                del_start    => $del_start,
                del_end      => $del_end,
                sample       => $sample,
                read_len     => $read_len,
                frag_n       => $wt_frag_n,
                insert_mean  => $insert_mean,
                insert_sd    => $insert_sd,
                error_rate   => $error_rate,
                r1_fh        => $R1,
                r2_fh        => $R2,
            );

            simulate_haplotype_reads(
                seq          => $del_seq,
                hap_name     => "DEL",
                chrom        => $chrom,
                gene_start   => $start,
                gene_end     => $end,
                del_start    => $del_start,
                del_end      => $del_end,
                sample       => $sample,
                read_len     => $read_len,
                frag_n       => $del_frag_n,
                insert_mean  => $insert_mean,
                insert_sd    => $insert_sd,
                error_rate   => $error_rate,
                r1_fh        => $R1,
                r2_fh        => $R2,
            );

            close $R1;
            close $R2;

            write_single_truth(
                file              => $truth_txt,
                sample            => $sample,
                gene              => $gene,
                transcript        => $transcript,
                chrom             => $chrom,
                strand            => $strand,
                gene_start        => $start,
                gene_end          => $end,
                gene_len          => $gene_len,
                exon_count        => $exon_count,
                bin_name          => $bin_name,
                min_size          => $min_size,
                max_size          => $max_size,
                rep               => $rep,
                del_start         => $del_start,
                del_end           => $del_end,
                del_size          => $del_size,
                wt_len            => $wt_len,
                del_hap_len       => $del_hap_len,
                depth             => $depth,
                wt_depth          => $wt_depth,
                del_depth         => $del_depth,
                read_len          => $read_len,
                insert_mean       => $insert_mean,
                insert_sd         => $insert_sd,
                error_rate        => $error_rate,
                wt_frag_n         => $wt_frag_n,
                del_frag_n        => $del_frag_n,
                seed              => $event_seed,
                r1_fastq          => $r1_fastq,
                r2_fastq          => $r2_fastq,
            );

            print $ALL join("\t",
                $sample,
                $gene,
                $transcript,
                $chrom,
                $strand,
                $start,
                $end,
                $gene_len,
                $exon_count,
                $bin_name,
                $min_size,
                $max_size,
                $rep,
                $del_start,
                $del_end,
                $del_size,
                $wt_len,
                $del_hap_len,
                $depth,
                $wt_depth,
                $del_depth,
                $read_len,
                $insert_mean,
                $insert_sd,
                $error_rate,
                $wt_frag_n,
                $del_frag_n,
                $event_seed,
                $r1_fastq,
                $r2_fastq
            ), "\n";

            print STDERR "[INFO] Generated $sample: $gene $chrom:$del_start-$del_end size=$del_size bp bin=$bin_name\n";
        }
    }
}

close $IN;
close $ALL;

print STDERR "[INFO] All truth file       : $all_truth\n";
print STDERR "[INFO] Genes processed      : $gene_count\n";
print STDERR "[INFO] Events generated     : $event_count\n";
print STDERR "[INFO] Events skipped       : $skip_count\n";
print STDERR "[INFO] Done.\n";

exit 0;


# ============================================================
# Functions
# ============================================================

sub clean_field {
    my ($x) = @_;
    $x = "" unless defined $x;
    $x =~ s/^\s+|\s+$//g;
    return $x;
}


sub sanitize_name {
    my ($x) = @_;
    $x =~ s/[^A-Za-z0-9_.-]/_/g;
    return $x;
}


sub random_int {
    my ($min, $max) = @_;

    die "[ERROR] Invalid random range: min=$min max=$max\n"
        if $min > $max;

    return int(rand($max - $min + 1)) + $min;
}


sub fetch_fasta_sequence {
    my ($samtools, $ref, $chrom, $start, $end) = @_;

    my $region = "$chrom:$start-$end";

    open my $FA, "-|", $samtools, "faidx", $ref, $region
        or die "[ERROR] Failed to run: $samtools faidx $ref $region\n";

    my $seq = "";
    while (my $line = <$FA>) {
        chomp $line;
        next if $line =~ /^>/;
        $line =~ s/\s+//g;
        $seq .= uc($line);
    }

    close $FA;

    die "[ERROR] Empty sequence extracted for region: $region\n"
        if $seq eq "";

    return $seq;
}


sub simulate_haplotype_reads {
    my %args = @_;

    my $seq         = $args{seq};
    my $hap_name    = $args{hap_name};
    my $chrom       = $args{chrom};
    my $gene_start  = $args{gene_start};
    my $gene_end    = $args{gene_end};
    my $del_start   = $args{del_start};
    my $del_end     = $args{del_end};
    my $sample      = $args{sample};
    my $read_len    = $args{read_len};
    my $frag_n      = $args{frag_n};
    my $insert_mean = $args{insert_mean};
    my $insert_sd   = $args{insert_sd};
    my $error_rate  = $args{error_rate};
    my $R1          = $args{r1_fh};
    my $R2          = $args{r2_fh};

    my $hap_len = length($seq);

    die "[ERROR] Haplotype length $hap_len < read length $read_len for $sample $hap_name\n"
        if $hap_len < $read_len;

    for my $i (1 .. $frag_n) {

        my $insert = random_insert_size($insert_mean, $insert_sd, $read_len, $hap_len);

        my $max_start0 = $hap_len - $insert;
        my $frag_start0 = int(rand($max_start0 + 1));

        my $fragment = substr($seq, $frag_start0, $insert);

        my $r1_seq = substr($fragment, 0, $read_len);
        my $r2_seq = substr($fragment, $insert - $read_len, $read_len);
        $r2_seq = revcomp($r2_seq);

        $r1_seq = introduce_errors($r1_seq, $error_rate);
        $r2_seq = introduce_errors($r2_seq, $error_rate);

        my $qual = "I" x $read_len;

        my $read_id = join(":",
            $sample,
            "hap=$hap_name",
            "frag=$i",
            "region=$chrom-$gene_start-$gene_end",
            "del=$del_start-$del_end",
            "hap_start0=$frag_start0",
            "insert=$insert"
        );

        print $R1 "\@$read_id/1\n";
        print $R1 "$r1_seq\n";
        print $R1 "+\n";
        print $R1 "$qual\n";

        print $R2 "\@$read_id/2\n";
        print $R2 "$r2_seq\n";
        print $R2 "+\n";
        print $R2 "$qual\n";
    }
}


sub random_insert_size {
    my ($mean, $sd, $min_insert, $max_insert) = @_;

    $min_insert = 1 if $min_insert < 1;

    die "[ERROR] max_insert $max_insert < min_insert $min_insert\n"
        if $max_insert < $min_insert;

    for my $try (1 .. 100) {
        my $x = int(rand_normal($mean, $sd) + 0.5);
        next if $x < $min_insert;
        next if $x > $max_insert;
        return $x;
    }

    my $fallback = int($mean + 0.5);
    $fallback = $min_insert if $fallback < $min_insert;
    $fallback = $max_insert if $fallback > $max_insert;

    return $fallback;
}


sub rand_normal {
    my ($mean, $sd) = @_;

    my $u1 = rand();
    my $u2 = rand();

    $u1 = 1e-10 if $u1 <= 0;

    my $z = sqrt(-2 * log($u1)) * cos(2 * 3.141592653589793 * $u2);

    return $mean + $sd * $z;
}


sub revcomp {
    my ($seq) = @_;
    $seq = reverse $seq;
    $seq =~ tr/ACGTNacgtn/TGCANtgcan/;
    return uc($seq);
}


sub introduce_errors {
    my ($seq, $error_rate) = @_;

    return $seq if $error_rate <= 0;

    my @bases = split //, uc($seq);
    my @dna = qw/A C G T/;

    for my $i (0 .. $#bases) {
        next if rand() >= $error_rate;

        my $b = $bases[$i];

        if ($b =~ /^[ACGT]$/) {
            my @alt = grep { $_ ne $b } @dna;
            $bases[$i] = $alt[int(rand(@alt))];
        }
    }

    return join("", @bases);
}


sub write_single_truth {
    my %x = @_;

    open my $T, ">", $x{file}
        or die "[ERROR] Cannot write truth file: $x{file}\n";

    print $T "Sample\t$x{sample}\n";
    print $T "Gene\t$x{gene}\n";
    print $T "Transcript\t$x{transcript}\n";
    print $T "Chrom\t$x{chrom}\n";
    print $T "Strand\t$x{strand}\n";
    print $T "Gene_Start\t$x{gene_start}\n";
    print $T "Gene_End\t$x{gene_end}\n";
    print $T "Gene_Length\t$x{gene_len}\n";
    print $T "ExonCount\t$x{exon_count}\n";
    print $T "Length_Bin\t$x{bin_name}\n";
    print $T "Bin_Min_Size\t$x{min_size}\n";
    print $T "Bin_Max_Size\t$x{max_size}\n";
    print $T "Replicate\t$x{rep}\n";
    print $T "Del_Start\t$x{del_start}\n";
    print $T "Del_End\t$x{del_end}\n";
    print $T "Del_Size_bp\t$x{del_size}\n";
    print $T "WT_Haplotype_Length\t$x{wt_len}\n";
    print $T "DEL_Haplotype_Length\t$x{del_hap_len}\n";
    print $T "Depth\t$x{depth}\n";
    print $T "WT_Depth\t$x{wt_depth}\n";
    print $T "DEL_Depth\t$x{del_depth}\n";
    print $T "Read_Length\t$x{read_len}\n";
    print $T "Insert_Mean\t$x{insert_mean}\n";
    print $T "Insert_SD\t$x{insert_sd}\n";
    print $T "Error_Rate\t$x{error_rate}\n";
    print $T "WT_Fragment_Number\t$x{wt_frag_n}\n";
    print $T "DEL_Fragment_Number\t$x{del_frag_n}\n";
    print $T "Seed\t$x{seed}\n";
    print $T "R1_FASTQ\t$x{r1_fastq}\n";
    print $T "R2_FASTQ\t$x{r2_fastq}\n";

    close $T;
}


sub usage {
    return <<"USAGE";

Usage:
  perl simulate_random_gene_wgs_del.pl \\
    --ref hg19.fa \\
    --gene-txt RefSeq_MANE_Select.gene.txt \\
    --outdir random_gene_del_sim \\
    --samtools samtools \\
    --replicates 3 \\
    --seed 20260527 \\
    --sample-prefix HCM \\
    --read-len 150 \\
    --depth 30 \\
    --insert-mean 350 \\
    --insert-sd 50 \\
    --error-rate 0.001

Required:
  --ref               Reference FASTA file. Must be indexed by samtools faidx.
  --gene-txt          Gene coordinate file.
                      Required columns:
                      Gene Transcript Chrom Start End Strand ExonCount

Optional:
  --gene-list         Optional gene list file. First column should be gene symbol.
                      If provided, only genes in this list will be simulated.

  --outdir            Output directory.
                      Default: random_gene_del_sim

  --samtools          Path to samtools.
                      Default: samtools

  --replicates        Number of random deletions per gene per size bin.
                      Default: 3

  --seed              Base random seed.
                      Default: 20260527

  --sample-prefix     Prefix for sample names.
                      Default: HCM

  --read-len          Read length.
                      Default: 150

  --depth             Total simulated depth.
                      WT haplotype uses depth / 2.
                      DEL haplotype uses depth / 2.
                      Default: 30

  --insert-mean       Mean insert size.
                      Default: 350

  --insert-sd         Insert size standard deviation.
                      Default: 50

  --error-rate        Per-base sequencing error rate.
                      Default: 0.001

  --max-genes         Debug option. Process only the first N genes.
                      Default: 0, meaning no limit.

Deletion size bins:
  50-100 bp
  100-200 bp
  200-500 bp
  500-1000 bp
  1000-2000 bp
  2000-5000 bp
  5000-10000 bp

Outputs:
  outdir/random_gene_del.all.truth.tsv
  outdir/<Sample>/<Sample>.R1.fastq
  outdir/<Sample>/<Sample>.R2.fastq
  outdir/<Sample>/<Sample>.truth.txt

USAGE
}

