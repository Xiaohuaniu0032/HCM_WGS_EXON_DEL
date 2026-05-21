#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);

# ============================================================
# Script: extract_sa_split_reads.pl
#
# Purpose:
#   Extract primary paired reads with SA tag from BAM and infer
#   deletion-like SV intervals from primary and SA aligned blocks.
#
# Required parameters:
#   --conf     Config file
#   --bam      Input BAM file
#   --out      Output TSV file
#
# Optional parameters:
#   --sample   Sample ID. If not provided, BAM basename without .bam is used.
#
# Region selection:
#   If ANALYZE_CORE_GENES_ONLY=1:
#       Only scan genes listed in HCM_CORE_GENE_LIST.
#       Gene coordinates are loaded from REFSEQ_MANE_SELECT_GENE_TXT.
#
#   If ANALYZE_CORE_GENES_ONLY=0:
#       Scan the whole BAM.
#
# Filtering rule:
#   Keep reads that are:
#     1. PAIRED
#     2. not UNMAPPED
#     3. not SECONDARY
#     4. not SUPPLEMENTARY
#     5. with SA tag
#
# SV interval definition:
#   Primary block = [Primary_Start, Primary_End]
#   SA block      = [SA_Start, SA_End]
#
#   For same-chromosome, non-overlapping blocks:
#     Left_Block  = block with smaller coordinate
#     Right_Block = block with larger coordinate
#
#     Left_Breakpoint  = Left_Block_End
#     Right_Breakpoint = Right_Block_Start
#
#     SV_Start  = Left_Block_End + 1
#     SV_End    = Right_Block_Start - 1
#     SV_Length = SV_End - SV_Start + 1
#
# Output:
#   Compact TSV file for downstream SA split-read clustering.
#   Primary_MAPQ and SA_MAPQ are retained for quality filtering.
# ============================================================

my $conf   = "";
my $bam    = "";
my $out    = "";
my $sample = "";
my $help   = 0;

GetOptions(
    "conf=s"   => \$conf,
    "bam=s"    => \$bam,
    "out=s"    => \$out,
    "sample=s" => \$sample,
    "help"     => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless $conf && $bam && $out;

die "[ERROR] Config file not found: $conf\n" unless -s $conf;
die "[ERROR] BAM file not found: $bam\n" unless -s $bam;

my $project_root = get_project_root($conf);
my %cfg = read_conf($conf);

my $samtools = get_conf_value(
    \%cfg,
    ["SAMTOOLS", "SAMTOOLS_BIN", "SAMTOOLS_PATH"],
    "samtools"
);

my $min_mapq = get_conf_value(
    \%cfg,
    ["SPLIT_READ_MIN_MAPQ", "MIN_SPLIT_READ_MAPQ", "MIN_MAPQ"],
    0
);

my $exclude_dup = get_conf_value(
    \%cfg,
    ["SPLIT_READ_EXCLUDE_DUP", "EXCLUDE_DUP", "EXCLUDE_DUPLICATE"],
    0
);

my $analyze_core_genes_only = get_conf_value(
    \%cfg,
    ["ANALYZE_CORE_GENES_ONLY"],
    0
);

if (!$sample) {
    $sample = basename($bam);
    $sample =~ s/\.bam$//;
}

my @scan_regions;

if ($analyze_core_genes_only) {
    my $core_gene_list = get_conf_value(
        \%cfg,
        ["HCM_CORE_GENE_LIST"],
        ""
    );

    my $gene_txt = get_conf_value(
        \%cfg,
        ["REFSEQ_MANE_SELECT_GENE_TXT", "REFSEQ_MANE_GENE_TXT", "GENE_TXT"],
        ""
    );

    die "[ERROR] ANALYZE_CORE_GENES_ONLY=1, but HCM_CORE_GENE_LIST is not configured\n"
        unless $core_gene_list;

    die "[ERROR] ANALYZE_CORE_GENES_ONLY=1, but REFSEQ_MANE_SELECT_GENE_TXT is not configured\n"
        unless $gene_txt;

    $core_gene_list = resolve_path($core_gene_list, $project_root);
    $gene_txt       = resolve_path($gene_txt,       $project_root);

    die "[ERROR] HCM_CORE_GENE_LIST file not found: $core_gene_list\n"
        unless -s $core_gene_list;

    die "[ERROR] REFSEQ_MANE_SELECT_GENE_TXT file not found: $gene_txt\n"
        unless -s $gene_txt;

    my %core_genes   = load_core_gene_list($core_gene_list);
    my %gene_regions = load_gene_regions($gene_txt);

    foreach my $gene (sort keys %core_genes) {
        unless (exists $gene_regions{$gene}) {
            warn "[WARN] Core gene not found in gene coordinate file: $gene\n";
            next;
        }

        foreach my $r (@{ $gene_regions{$gene} }) {
            push @scan_regions, $r;
        }
    }

    die "[ERROR] No scan regions generated from HCM_CORE_GENE_LIST and REFSEQ_MANE_SELECT_GENE_TXT\n"
        unless @scan_regions;
}
else {
    push @scan_regions, {
        gene       => "ALL",
        transcript => "NA",
        chrom      => "ALL",
        start      => "NA",
        end        => "NA",
        region     => "",
        name       => "Whole_BAM",
    };
}

open my $out_fh, ">", $out
    or die "[ERROR] Cannot write output file: $out\n";

print $out_fh join("\t", qw(
    Sample
    Scan_Mode
    Scan_Gene
    Scan_Region
    Read_Name
    Flag
    Primary_Chrom
    Primary_Start
    Primary_End
    Primary_CIGAR
    Primary_MAPQ
    Primary_Strand
    SA_Chrom
    SA_Start
    SA_End
    SA_CIGAR
    SA_MAPQ
    SA_Strand
    Left_Breakpoint
    Right_Breakpoint
    SV_Chrom
    SV_Start
    SV_End
    SV_Length
    Junction_Type
    SA_Tag
)) . "\n";

my $scan_mode = $analyze_core_genes_only ? "Core_Genes" : "Whole_BAM";

my $total_records     = 0;
my $qualified_reads   = 0;
my $qualified_events  = 0;
my $skipped_no_sa     = 0;
my $skipped_flag      = 0;
my $skipped_mapq      = 0;
my $skipped_bad_cigar = 0;

foreach my $scan_region (@scan_regions) {
    my $region_string = $scan_region->{region};

    my $cmd = build_samtools_view_cmd(
        samtools => $samtools,
        bam      => $bam,
        region   => $region_string,
    );

    open my $pipe, "$cmd |"
        or die "[ERROR] Failed to run command: $cmd\n";

    while (my $line = <$pipe>) {
        chomp $line;
        next if $line =~ /^\s*$/;

        $total_records++;

        my @f = split /\t/, $line;
        next unless @f >= 11;

        my ($qname, $flag, $rname, $pos, $mapq, $cigar) = @f[0..5];

        unless (defined $flag && $flag =~ /^\d+$/) {
            $skipped_flag++;
            next;
        }

        unless (
            is_paired($flag)
            && !is_unmapped($flag)
            && !is_secondary($flag)
            && !is_supplementary($flag)
        ) {
            $skipped_flag++;
            next;
        }

        if ($exclude_dup && is_duplicate($flag)) {
            $skipped_flag++;
            next;
        }

        unless (defined $mapq && $mapq =~ /^\d+$/ && $mapq >= $min_mapq) {
            $skipped_mapq++;
            next;
        }

        unless ($cigar && $cigar ne "*") {
            $skipped_bad_cigar++;
            next;
        }

        my ($has_sa, $sa_tag) = parse_sa_tag(@f);
        unless ($has_sa) {
            $skipped_no_sa++;
            next;
        }

        my $primary_ref_len = cigar_ref_length($cigar);
        unless ($primary_ref_len > 0) {
            $skipped_bad_cigar++;
            next;
        }

        my $primary_start  = $pos;
        my $primary_end    = $pos + $primary_ref_len - 1;
        my $primary_strand = is_reverse($flag) ? "-" : "+";

        my @sa_entries = parse_sa_entries($sa_tag);
        next unless @sa_entries;

        $qualified_reads++;

        foreach my $sa (@sa_entries) {
            next unless defined $sa->{chr} && $sa->{chr} ne "";
            next unless defined $sa->{pos} && $sa->{pos} =~ /^\d+$/;
            next unless defined $sa->{cigar} && $sa->{cigar} ne "*";

            my $sa_ref_len = cigar_ref_length($sa->{cigar});
            next unless $sa_ref_len > 0;

            my $sa_start = $sa->{pos};
            my $sa_end   = $sa->{pos} + $sa_ref_len - 1;

            my %sv = infer_sv_interval(
                primary_chr   => $rname,
                primary_start => $primary_start,
                primary_end   => $primary_end,
                sa_chr        => $sa->{chr},
                sa_start      => $sa_start,
                sa_end        => $sa_end,
            );

            print $out_fh join("\t",
                $sample,
                $scan_mode,
                $scan_region->{gene},
                $scan_region->{region} || "Whole_BAM",

                $qname,
                $flag,

                $rname,
                $primary_start,
                $primary_end,
                $cigar,
                $mapq,
                $primary_strand,

                $sa->{chr},
                $sa_start,
                $sa_end,
                $sa->{cigar},
                $sa->{mapq},
                $sa->{strand},

                $sv{left_breakpoint},
                $sv{right_breakpoint},
                $sv{sv_chr},
                $sv{sv_start},
                $sv{sv_end},
                $sv{sv_length},
                $sv{junction_type},

                $sa_tag,
            ) . "\n";

            $qualified_events++;
        }
    }

    close $pipe;
}

close $out_fh;

print STDERR "[INFO] SA split-read extraction finished\n";
print STDERR "[INFO] Config                  : $conf\n";
print STDERR "[INFO] Project root            : $project_root\n";
print STDERR "[INFO] BAM                     : $bam\n";
print STDERR "[INFO] Output                  : $out\n";
print STDERR "[INFO] Sample                  : $sample\n";
print STDERR "[INFO] samtools                : $samtools\n";
print STDERR "[INFO] ANALYZE_CORE_GENES_ONLY : $analyze_core_genes_only\n";
print STDERR "[INFO] Scan mode               : $scan_mode\n";
print STDERR "[INFO] Scan region number      : " . scalar(@scan_regions) . "\n";
print STDERR "[INFO] Minimum primary MAPQ    : $min_mapq\n";
print STDERR "[INFO] Exclude duplicates      : $exclude_dup\n";
print STDERR "[INFO] Total SAM records       : $total_records\n";
print STDERR "[INFO] Qualified reads         : $qualified_reads\n";
print STDERR "[INFO] Qualified SA events     : $qualified_events\n";
print STDERR "[INFO] Skipped by flag         : $skipped_flag\n";
print STDERR "[INFO] Skipped by MAPQ         : $skipped_mapq\n";
print STDERR "[INFO] Skipped without SA      : $skipped_no_sa\n";
print STDERR "[INFO] Skipped bad CIGAR       : $skipped_bad_cigar\n";

exit 0;


sub usage {
    return <<"USAGE";
Usage:
  perl extract_sa_split_reads.pl \\
      --conf   hcm_exondel.example.conf \\
      --bam    input.bam \\
      --out    output.sa_split_reads.tsv \\
      --sample SAMPLE_ID

Required:
  --conf STR      Config file
  --bam  STR      Input BAM file
  --out  STR      Output TSV file

Optional:
  --sample STR    Sample ID. Default: BAM basename without .bam

Config keys used:
  SAMTOOLS / SAMTOOLS_BIN / SAMTOOLS_PATH

  ANALYZE_CORE_GENES_ONLY
      1: scan only genes in HCM_CORE_GENE_LIST
      0: scan whole BAM

  HCM_CORE_GENE_LIST
      Required when ANALYZE_CORE_GENES_ONLY=1

  REFSEQ_MANE_SELECT_GENE_TXT
      Required when ANALYZE_CORE_GENES_ONLY=1

  SPLIT_READ_MIN_MAPQ / MIN_SPLIT_READ_MAPQ / MIN_MAPQ
      Minimum MAPQ for primary alignment.

  SPLIT_READ_EXCLUDE_DUP / EXCLUDE_DUP / EXCLUDE_DUPLICATE

Example:
  perl bin/extract_sa_split_reads.pl \\
      --conf conf/hcm_exondel.example.conf \\
      --bam test/test_results/25B09089386.final.merge/05.gene_bam/FHOD3.bam \\
      --out test/test_results/25B09089386.final.merge/05.gene_bam/FHOD3.sa_split_reads.tsv \\
      --sample 25B09089386

USAGE
}


sub get_project_root {
    my ($conf_file) = @_;

    my $abs_conf = abs_path($conf_file);
    my $dir = dirname($abs_conf);

    if (basename($dir) eq "conf") {
        return dirname($dir);
    }

    return $dir;
}


sub resolve_path {
    my ($path, $project_root) = @_;

    return $path if $path =~ m{^/};

    my $full = "$project_root/$path";
    $full =~ s{//+}{/}g;

    return $full;
}


sub read_conf {
    my ($file) = @_;

    my %cfg;

    open my $fh, "<", $file
        or die "[ERROR] Cannot open config file: $file\n";

    while (my $line = <$fh>) {
        chomp $line;

        $line =~ s/\r$//;
        $line =~ s/^\s+//;
        $line =~ s/\s+$//;

        next if $line eq "";
        next if $line =~ /^#/;

        $line =~ s/\s+#.*$//;

        next unless $line =~ /^([A-Za-z0-9_.-]+)\s*=\s*(.*)$/;

        my $key = $1;
        my $val = $2;

        $val =~ s/^\s+//;
        $val =~ s/\s+$//;
        $val =~ s/^["']//;
        $val =~ s/["']$//;

        $cfg{$key} = $val;
    }

    close $fh;

    return %cfg;
}


sub get_conf_value {
    my ($cfg_ref, $keys_ref, $default) = @_;

    foreach my $key (@$keys_ref) {
        if (exists $cfg_ref->{$key} && defined $cfg_ref->{$key} && $cfg_ref->{$key} ne "") {
            return $cfg_ref->{$key};
        }
    }

    return $default;
}


sub load_core_gene_list {
    my ($file) = @_;

    my %genes;

    open my $fh, "<", $file
        or die "[ERROR] Cannot open HCM core gene list: $file\n";

    while (my $line = <$fh>) {
        chomp $line;

        $line =~ s/\r$//;
        $line =~ s/^\s+//;
        $line =~ s/\s+$//;

        next if $line eq "";
        next if $line =~ /^#/;

        my @f = split /\t/, $line;
        my $gene = $f[0];

        $gene =~ s/^\s+//;
        $gene =~ s/\s+$//;

        next if $gene eq "";
        next if uc($gene) eq "GENE";
        next if uc($gene) eq "GENE_SYMBOL";

        $genes{$gene} = 1;
    }

    close $fh;

    return %genes;
}


sub load_gene_regions {
    my ($file) = @_;

    my %regions;

    open my $fh, "<", $file
        or die "[ERROR] Cannot open gene coordinate file: $file\n";

    my $header = <$fh>;
    die "[ERROR] Empty gene coordinate file: $file\n" unless defined $header;

    chomp $header;
    $header =~ s/\r$//;

    my @header = split /\t/, $header;
    my %idx;

    for (my $i = 0; $i < @header; $i++) {
        my $col = $header[$i];
        $col =~ s/^\s+//;
        $col =~ s/\s+$//;
        $idx{$col} = $i;
    }

    my @gene_keys  = qw(Gene gene GeneSymbol gene_symbol Symbol symbol);
    my @tx_keys    = qw(Transcript transcript Transcript_ID transcript_id);
    my @chr_keys   = qw(Chrom chrom Chr chr Chromosome chromosome);
    my @start_keys = qw(Start start Gene_Start gene_start);
    my @end_keys   = qw(End end Gene_End gene_end);

    my $gene_idx  = find_col_index(\%idx, \@gene_keys);
    my $tx_idx    = find_col_index(\%idx, \@tx_keys);
    my $chr_idx   = find_col_index(\%idx, \@chr_keys);
    my $start_idx = find_col_index(\%idx, \@start_keys);
    my $end_idx   = find_col_index(\%idx, \@end_keys);

    my $has_header = (
        defined $gene_idx
        && defined $chr_idx
        && defined $start_idx
        && defined $end_idx
    ) ? 1 : 0;

    if (!$has_header) {
        seek($fh, 0, 0);
    }

    while (my $line = <$fh>) {
        chomp $line;

        $line =~ s/\r$//;
        next if $line =~ /^\s*$/;
        next if $line =~ /^#/;

        my @f = split /\t/, $line;

        my ($gene, $tx, $chr, $start, $end);

        if ($has_header) {
            $gene  = safe_field(\@f, $gene_idx);
            $tx    = defined $tx_idx ? safe_field(\@f, $tx_idx) : "NA";
            $chr   = safe_field(\@f, $chr_idx);
            $start = safe_field(\@f, $start_idx);
            $end   = safe_field(\@f, $end_idx);
        }
        else {
            # Common format 1:
            # Gene Transcript Chrom Start End
            if (@f >= 5 && $f[2] =~ /^chr|^[0-9XYM]+$/i && $f[3] =~ /^\d+$/ && $f[4] =~ /^\d+$/) {
                ($gene, $tx, $chr, $start, $end) = @f[0..4];
            }
            # Common format 2:
            # Chrom Start End Gene Transcript
            elsif (@f >= 5 && $f[1] =~ /^\d+$/ && $f[2] =~ /^\d+$/) {
                ($chr, $start, $end, $gene, $tx) = @f[0..4];
            }
            else {
                next;
            }
        }

        next unless defined $gene && defined $chr && defined $start && defined $end;
        next unless $start =~ /^\d+$/ && $end =~ /^\d+$/;

        if ($start > $end) {
            my $tmp = $start;
            $start = $end;
            $end = $tmp;
        }

        my $region = "$chr:$start-$end";
        my $name   = join("|", $gene, $tx || "NA", $region);

        push @{ $regions{$gene} }, {
            gene       => $gene,
            transcript => $tx || "NA",
            chrom      => $chr,
            start      => $start,
            end        => $end,
            region     => $region,
            name       => $name,
        };
    }

    close $fh;

    return %regions;
}


sub find_col_index {
    my ($idx_ref, $keys_ref) = @_;

    foreach my $key (@$keys_ref) {
        return $idx_ref->{$key} if exists $idx_ref->{$key};
    }

    return undef;
}


sub safe_field {
    my ($arr_ref, $idx) = @_;

    return "NA" unless defined $idx;
    return "NA" unless defined $arr_ref->[$idx];

    my $v = $arr_ref->[$idx];
    $v =~ s/^\s+//;
    $v =~ s/\s+$//;

    return $v eq "" ? "NA" : $v;
}


sub build_samtools_view_cmd {
    my %args = @_;

    my $samtools = $args{samtools};
    my $bam      = $args{bam};
    my $region   = $args{region};

    my @cmd = (
        shell_quote($samtools),
        "view",
        shell_quote($bam),
    );

    if (defined $region && $region ne "") {
        push @cmd, shell_quote($region);
    }

    return join(" ", @cmd);
}


sub parse_sa_tag {
    my @fields = @_;

    for (my $i = 11; $i < @fields; $i++) {
        if ($fields[$i] =~ /^SA:Z:(.+)$/) {
            return (1, $1);
        }
    }

    return (0, "NA");
}


sub parse_sa_entries {
    my ($sa_tag) = @_;

    my @entries;
    return @entries if !$sa_tag || $sa_tag eq "NA";

    my @items = split /;/, $sa_tag;

    foreach my $item (@items) {
        next if $item =~ /^\s*$/;

        my @f = split /,/, $item;

        # SA format:
        # rname,pos,strand,CIGAR,mapQ,NM;
        next unless @f >= 6;

        my ($chr, $pos, $strand, $cigar, $mapq, $nm) = @f[0..5];

        $mapq = "NA" unless defined $mapq && $mapq ne "";

        push @entries, {
            chr    => $chr,
            pos    => $pos,
            strand => $strand,
            cigar  => $cigar,
            mapq   => $mapq,
            nm     => $nm,
        };
    }

    return @entries;
}


sub cigar_ref_length {
    my ($cigar) = @_;

    return 0 if !$cigar || $cigar eq "*";

    my $len = 0;

    while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
        my ($n, $op) = ($1, $2);

        # Operations consuming reference:
        # M, D, N, =, X
        if ($op =~ /^[MDN=X]$/) {
            $len += $n;
        }
    }

    return $len;
}


sub infer_sv_interval {
    my %args = @_;

    my $primary_chr   = $args{primary_chr};
    my $primary_start = $args{primary_start};
    my $primary_end   = $args{primary_end};

    my $sa_chr        = $args{sa_chr};
    my $sa_start      = $args{sa_start};
    my $sa_end        = $args{sa_end};

    my %ret = (
        left_breakpoint   => "NA",
        right_breakpoint  => "NA",
        sv_chr            => "NA",
        sv_start          => "NA",
        sv_end            => "NA",
        sv_length         => "NA",
        junction_type     => "NA",
    );

    if ($primary_chr ne $sa_chr) {
        $ret{junction_type} = "Interchromosomal";
        return %ret;
    }

    my ($left_start, $left_end, $right_start, $right_end);

    if ($primary_start <= $sa_start) {
        ($left_start,  $left_end)  = ($primary_start, $primary_end);
        ($right_start, $right_end) = ($sa_start,      $sa_end);
    }
    else {
        ($left_start,  $left_end)  = ($sa_start,      $sa_end);
        ($right_start, $right_end) = ($primary_start, $primary_end);
    }

    $ret{left_breakpoint}  = $left_end;
    $ret{right_breakpoint} = $right_start;

    if ($left_end < $right_start - 1) {
        my $sv_start  = $left_end + 1;
        my $sv_end    = $right_start - 1;
        my $sv_length = $sv_end - $sv_start + 1;

        $ret{sv_chr}        = $primary_chr;
        $ret{sv_start}      = $sv_start;
        $ret{sv_end}        = $sv_end;
        $ret{sv_length}     = $sv_length;
        $ret{junction_type} = "Deletion_Gap";
    }
    elsif ($left_end == $right_start - 1) {
        $ret{sv_chr}        = $primary_chr;
        $ret{sv_start}      = "NA";
        $ret{sv_end}        = "NA";
        $ret{sv_length}     = 0;
        $ret{junction_type} = "Adjacent_No_Gap";
    }
    else {
        $ret{sv_chr}        = $primary_chr;
        $ret{sv_start}      = "NA";
        $ret{sv_end}        = "NA";
        $ret{sv_length}     = "NA";
        $ret{junction_type} = "Overlap_or_Complex";
    }

    return %ret;
}


sub is_paired {
    my ($flag) = @_;
    return ($flag & 0x1) ? 1 : 0;
}


sub is_unmapped {
    my ($flag) = @_;
    return ($flag & 0x4) ? 1 : 0;
}


sub is_reverse {
    my ($flag) = @_;
    return ($flag & 0x10) ? 1 : 0;
}


sub is_secondary {
    my ($flag) = @_;
    return ($flag & 0x100) ? 1 : 0;
}


sub is_duplicate {
    my ($flag) = @_;
    return ($flag & 0x400) ? 1 : 0;
}


sub is_supplementary {
    my ($flag) = @_;
    return ($flag & 0x800) ? 1 : 0;
}


sub shell_quote {
    my ($s) = @_;

    die "[ERROR] Undefined shell argument\n" unless defined $s;

    $s =~ s/'/'"'"'/g;
    return "'$s'";
}

