#!/usr/bin/env perl
# -*- coding: utf-8 -*-

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(basename);
use File::Path qw(make_path);
use Cwd qw(abs_path);
use FindBin qw($Bin);

# ============================================================
# HCMExonDel main pipeline
#
# Design:
#   1. Only BAM list mode is supported.
#   2. One BAM file path per line.
#   3. Sample name is derived from BAM basename.
#   4. All interval files use tab-delimited TXT format.
#   5. Coordinates in interval TXT files are 1-based closed.
#   6. No BED conversion is performed in this main script.
#
# Required interval files in config:
#   REFSEQ_MANE_SELECT_EXON_TXT
#   REFSEQ_MANE_SELECT_GENE_TXT
#
# Workflow:
#   1. gene_mean_depth
#   2. split_reads
#   3. discordant_reads
#   4. merge_evidence
#   5. annotate_candidates
#
# Note:
#   check_bam.pl is not called at this stage.
# ============================================================

my $config;
my $bam_list;
my $outdir = "results";
my $shell;
my $force = 0;
my $help  = 0;

GetOptions(
    "config=s"   => \$config,
    "bam-list=s" => \$bam_list,
    "outdir=s"   => \$outdir,
    "shell=s"    => \$shell,
    "force"      => \$force,
    "help"       => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless $config && $bam_list;

$config = abs_path($config);
die "[ERROR] Config file not found: $config\n" unless $config && -s $config;

$bam_list = abs_path($bam_list);
die "[ERROR] BAM list file not found: $bam_list\n" unless $bam_list && -s $bam_list;

# -----------------------------
# Read config
# -----------------------------
my %CONF = read_config($config);

my $perl     = get_conf(\%CONF, "PERL", "perl");
my $threads  = get_conf(\%CONF, "THREADS", 4);
my $keep_tmp = get_conf(\%CONF, "KEEP_TMP", 0);

my $exon_txt = get_conf(\%CONF, "REFSEQ_MANE_SELECT_EXON_TXT", "");
my $gene_txt = get_conf(\%CONF, "REFSEQ_MANE_SELECT_GENE_TXT", "");

die "[ERROR] REFSEQ_MANE_SELECT_EXON_TXT is required in config\n" unless $exon_txt;
die "[ERROR] REFSEQ_MANE_SELECT_GENE_TXT is required in config\n" unless $gene_txt;

$exon_txt = abs_path($exon_txt);
$gene_txt = abs_path($gene_txt);

die "[ERROR] Exon TXT file not found: $exon_txt\n" unless $exon_txt && -s $exon_txt;
die "[ERROR] Gene TXT file not found: $gene_txt\n" unless $gene_txt && -s $gene_txt;

check_exon_txt_format($exon_txt);
check_gene_txt_format($gene_txt);

# -----------------------------
# Prepare output directory
# -----------------------------
make_path($outdir) unless -d $outdir;
$outdir = abs_path($outdir);

$shell ||= "$outdir/run_hcm_exondel.sh";

if (-e $shell && !$force) {
    die "[ERROR] Shell file already exists: $shell\n"
      . "        Please use --force to overwrite it or specify another --shell\n";
}

# -----------------------------
# Check dependent scripts
# -----------------------------
my $script_depth      = "$Bin/bin/run_gene_mean_depth.pl";
my $script_split      = "$Bin/bin/run_split_reads.pl";
my $script_discordant = "$Bin/bin/run_discordant_reads.pl";
my $script_merge      = "$Bin/bin/merge_evidence.pl";
my $script_annotate   = "$Bin/bin/annotate_candidates.pl";

check_script($script_depth);
check_script($script_split);
check_script($script_discordant);
check_script($script_merge);
check_script($script_annotate);

# -----------------------------
# Read BAM list
# -----------------------------
my @samples = read_bam_list($bam_list);
die "[ERROR] No valid BAM found in BAM list: $bam_list\n" unless @samples;

# -----------------------------
# Print summary
# -----------------------------
print "\n";
print "============================================================\n";
print " HCMExonDel shell generation\n";
print "============================================================\n";
print "Config file  : $config\n";
print "BAM list     : $bam_list\n";
print "Exon TXT     : $exon_txt\n";
print "Gene TXT     : $gene_txt\n";
print "Outdir       : $outdir\n";
print "Threads      : $threads\n";
print "Keep tmp     : $keep_tmp\n";
print "Samples      : " . scalar(@samples) . "\n";
print "Shell        : $shell\n";
print "============================================================\n\n";

# -----------------------------
# Generate shell
# -----------------------------
open my $SH, ">", $shell or die "[ERROR] Cannot write shell file: $shell\n";

print_shell_header(
    fh       => $SH,
    config   => $config,
    bam_list => $bam_list,
    exon_txt => $exon_txt,
    gene_txt => $gene_txt,
    outdir   => $outdir,
    threads  => $threads,
);

foreach my $item (@samples) {
    my $sample = $item->{sample};
    my $bam    = $item->{bam};

    my $sample_dir     = "$outdir/$sample";
    my $log_dir        = "$sample_dir/00.log";
    my $depth_dir      = "$sample_dir/01.depth";
    my $split_dir      = "$sample_dir/02.split_reads";
    my $discordant_dir = "$sample_dir/03.discordant_reads";
    my $candidate_dir  = "$sample_dir/04.candidates";
    my $report_dir     = "$sample_dir/05.report";
    my $tmp_dir        = "$sample_dir/tmp";

    if (-d $sample_dir && !$force) {
        die "[ERROR] Sample output directory already exists: $sample_dir\n"
          . "        Please use --force to allow existing sample directories\n";
    }

    make_path($log_dir);
    make_path($depth_dir);
    make_path($split_dir);
    make_path($discordant_dir);
    make_path($candidate_dir);
    make_path($report_dir);
    make_path($tmp_dir);

    my $depth_out      = "$depth_dir/$sample.depth_candidates.tsv";
    my $split_out      = "$split_dir/$sample.split_reads.tsv";
    my $discordant_out = "$discordant_dir/$sample.discordant_reads.tsv";
    my $merged_out     = "$candidate_dir/$sample.merged_candidates.tsv";
    my $annotated_out  = "$report_dir/$sample.annotated_candidates.tsv";

    print $SH "\n";
    print $SH "############################################################\n";
    print $SH "# Sample: $sample\n";
    print $SH "############################################################\n";
    print $SH "echo \"[INFO] Start sample: $sample \$(date)\"\n";
    print $SH "mkdir -p "
        . shell_quote($log_dir) . " "
        . shell_quote($depth_dir) . " "
        . shell_quote($split_dir) . " "
        . shell_quote($discordant_dir) . " "
        . shell_quote($candidate_dir) . " "
        . shell_quote($report_dir) . " "
        . shell_quote($tmp_dir) . "\n\n";

    print $SH "echo \"[INFO] Step 1: gene_mean_depth for $sample \$(date)\"\n";
    print $SH join(" ",
        shell_quote($perl),
        shell_quote($script_depth),
        "--config", shell_quote($config),
        "--bam",    shell_quote($bam),
        "--sample", shell_quote($sample),
        "--out",    shell_quote($depth_out),
        "--threads", $threads
    ) . " > " . shell_quote("$log_dir/01.gene_mean_depth.log") . " 2>&1\n\n";

    print $SH "echo \"[INFO] Step 2: split_reads for $sample \$(date)\"\n";
    print $SH join(" ",
        shell_quote($perl),
        shell_quote($script_split),
        "--config", shell_quote($config),
        "--bam",    shell_quote($bam),
        "--sample", shell_quote($sample),
        "--out",    shell_quote($split_out),
        "--threads", $threads
    ) . " > " . shell_quote("$log_dir/02.split_reads.log") . " 2>&1\n\n";

    print $SH "echo \"[INFO] Step 3: discordant_reads for $sample \$(date)\"\n";
    print $SH join(" ",
        shell_quote($perl),
        shell_quote($script_discordant),
        "--config", shell_quote($config),
        "--bam",    shell_quote($bam),
        "--sample", shell_quote($sample),
        "--out",    shell_quote($discordant_out),
        "--threads", $threads
    ) . " > " . shell_quote("$log_dir/03.discordant_reads.log") . " 2>&1\n\n";

    print $SH "echo \"[INFO] Step 4: merge_evidence for $sample \$(date)\"\n";
    print $SH join(" ",
        shell_quote($perl),
        shell_quote($script_merge),
        "--config",     shell_quote($config),
        "--sample",     shell_quote($sample),
        "--depth",      shell_quote($depth_out),
        "--split",      shell_quote($split_out),
        "--discordant", shell_quote($discordant_out),
        "--out",        shell_quote($merged_out)
    ) . " > " . shell_quote("$log_dir/04.merge_evidence.log") . " 2>&1\n\n";

    print $SH "echo \"[INFO] Step 5: annotate_candidates for $sample \$(date)\"\n";
    print $SH join(" ",
        shell_quote($perl),
        shell_quote($script_annotate),
        "--config", shell_quote($config),
        "--input",  shell_quote($merged_out),
        "--out",    shell_quote($annotated_out)
    ) . " > " . shell_quote("$log_dir/05.annotate_candidates.log") . " 2>&1\n\n";

    if (!$keep_tmp) {
        print $SH "rm -rf " . shell_quote($tmp_dir) . "\n\n";
    }

    print $SH "echo \"[INFO] Finished sample: $sample \$(date)\"\n";
    print $SH "echo \"[INFO] Annotated candidates: $annotated_out\"\n";
}

print $SH "\n";
print $SH "echo \"[INFO] All samples finished \$(date)\"\n";

close $SH;

chmod 0755, $shell or die "[ERROR] Cannot chmod shell file: $shell\n";

print "[INFO] Shell generated successfully\n";
print "[INFO] Shell file: $shell\n";
print "[INFO] Run command:\n";
print "       bash $shell\n\n";

exit 0;

# ============================================================
# Subroutines
# ============================================================

sub read_config {
    my ($file) = @_;
    my %conf;

    open my $FH, "<", $file or die "[ERROR] Cannot open config file: $file\n";

    while (my $line = <$FH>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        if ($line =~ /^\s*([^=\s]+)\s*=\s*(.*?)\s*$/) {
            my $key = $1;
            my $val = $2;

            $val =~ s/\s+#.*$//;
            $val =~ s/^\s+|\s+$//g;
            $val =~ s/^['"]//;
            $val =~ s/['"]$//;

            $conf{$key} = $val;
        }
    }

    close $FH;
    return %conf;
}

sub get_conf {
    my ($conf_ref, $key, $default) = @_;

    if (exists $conf_ref->{$key} && defined $conf_ref->{$key} && $conf_ref->{$key} ne "") {
        return $conf_ref->{$key};
    }

    return $default;
}

sub read_bam_list {
    my ($file) = @_;

    my @list;
    my %seen_sample;

    open my $FH, "<", $file or die "[ERROR] Cannot open BAM list: $file\n";

    my $line_no = 0;

    while (my $line = <$FH>) {
        chomp $line;
        $line =~ s/\r$//;
        $line_no++;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        my @fields = split /\s+/, $line;

        die "[ERROR] Invalid BAM list format at line $line_no: $line\n"
          . "        Current version requires one BAM file per line only.\n"
            unless @fields == 1;

        my $bam = $fields[0];
        my $abs_bam = abs_path($bam);

        die "[ERROR] BAM file not found at line $line_no: $bam\n"
            unless $abs_bam && -s $abs_bam;

        check_bam_index($abs_bam);

        my $sample = bam_to_sample_name($abs_bam);

        die "[ERROR] Empty sample name derived from BAM at line $line_no: $bam\n"
            unless $sample;

        if (exists $seen_sample{$sample}) {
            die "[ERROR] Duplicated sample name detected: $sample\n"
              . "        Please make BAM basenames unique.\n";
        }

        $seen_sample{$sample} = 1;

        push @list, {
            sample => $sample,
            bam    => $abs_bam,
        };
    }

    close $FH;
    return @list;
}

sub bam_to_sample_name {
    my ($bam) = @_;

    my $name = basename($bam);

    $name =~ s/\.sorted\.markdup\.bam$//i;
    $name =~ s/\.sorted\.dedup\.bam$//i;
    $name =~ s/\.markdup\.bam$//i;
    $name =~ s/\.dedup\.bam$//i;
    $name =~ s/\.recal\.bam$//i;
    $name =~ s/\.realn\.bam$//i;
    $name =~ s/\.sorted\.bam$//i;
    $name =~ s/\.bam$//i;

    $name =~ s/[^\w.\-]+/_/g;

    return $name;
}

sub check_bam_index {
    my ($bam) = @_;

    my $bai1 = "$bam.bai";
    my $bai2 = $bam;
    $bai2 =~ s/\.bam$/.bai/i;

    die "[ERROR] BAM index not found for: $bam\n"
      . "        Expected: $bai1 or $bai2\n"
        unless (-s $bai1 || -s $bai2);

    return 1;
}

sub check_script {
    my ($script) = @_;

    die "[ERROR] Required script not found: $script\n" unless -s $script;

    return 1;
}

sub check_exon_txt_format {
    my ($file) = @_;

    open my $FH, "<", $file or die "[ERROR] Cannot open exon TXT: $file\n";

    my $header = <$FH>;
    close $FH;

    chomp $header;
    $header =~ s/\r$//;

    my @h = split /\t/, $header;
    my %h = map { $_ => 1 } @h;

    foreach my $col (qw(Gene Transcript Exon Chrom Start End Strand)) {
        die "[ERROR] Exon TXT missing required column '$col': $file\n"
            unless exists $h{$col};
    }

    return 1;
}

sub check_gene_txt_format {
    my ($file) = @_;

    open my $FH, "<", $file or die "[ERROR] Cannot open gene TXT: $file\n";

    my $header = <$FH>;
    close $FH;

    chomp $header;
    $header =~ s/\r$//;

    my @h = split /\t/, $header;
    my %h = map { $_ => 1 } @h;

    foreach my $col (qw(Gene Transcript Chrom Start End Strand ExonCount)) {
        die "[ERROR] Gene TXT missing required column '$col': $file\n"
            unless exists $h{$col};
    }

    return 1;
}

sub print_shell_header {
    my %args = @_;

    my $FH       = $args{fh};
    my $config   = $args{config};
    my $bam_list = $args{bam_list};
    my $exon_txt = $args{exon_txt};
    my $gene_txt = $args{gene_txt};
    my $outdir   = $args{outdir};
    my $threads  = $args{threads};

    print $FH "#!/usr/bin/env bash\n";
    print $FH "set -euo pipefail\n\n";
    print $FH "############################################################\n";
    print $FH "# Auto-generated by hcm_exondel.pl\n";
    print $FH "# Config   : $config\n";
    print $FH "# BAM list : $bam_list\n";
    print $FH "# Exon TXT : $exon_txt\n";
    print $FH "# Gene TXT : $gene_txt\n";
    print $FH "# Outdir   : $outdir\n";
    print $FH "# Threads  : $threads\n";
    print $FH "############################################################\n\n";
    print $FH "echo \"[INFO] HCMExonDel workflow started \$(date)\"\n";
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
  perl hcm_exondel.pl \\
      --config conf/hcm_exondel.example.conf \\
      --bam-list example/input_bam.list \\
      --outdir results \\
      --force

Required arguments:
  --config       Configuration file.
  --bam-list     BAM list file. One BAM file path per line.

Optional arguments:
  --outdir       Output directory. [default: results]
  --shell        Output shell file. [default: <outdir>/run_hcm_exondel.sh]
  --force        Overwrite existing shell and allow existing sample directories.
  --help         Show this help message.

Required interval TXT files in config:
  REFSEQ_MANE_SELECT_EXON_TXT
  REFSEQ_MANE_SELECT_GENE_TXT

Expected exon TXT format:
  Gene    Transcript    Exon    Chrom    Start    End    Strand

Expected gene TXT format:
  Gene    Transcript    Chrom    Start    End    Strand    ExonCount

Coordinate:
  All TXT interval files are 1-based closed intervals.

BAM list format:
  /path/to/sample1.sorted.bam
  /path/to/sample2.markdup.bam
  /path/to/sample3.bam

Generated workflow:
  gene_mean_depth
  split_reads
  discordant_reads
  merge_evidence
  annotate_candidates

Output structure:
  results/
  ├── run_hcm_exondel.sh
  └── SAMPLE001/
      ├── 00.log/
      ├── 01.depth/
      ├── 02.split_reads/
      ├── 03.discordant_reads/
      ├── 04.candidates/
      ├── 05.report/
      └── tmp/

USAGE
}

