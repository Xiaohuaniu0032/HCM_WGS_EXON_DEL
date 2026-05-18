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
# Function:
#   Generate one shell script for each BAM sample.
#
# Input:
#   1. Config file
#   2. BAM list, one BAM file per line
#
# Output:
#   <outdir>/sample_shell.list
#   <outdir>/qsub_command.list
#   <outdir>/<sample>/<sample>.run.sh
#
# Workflow:
#   1. gene_mean_depth
#   2. split_reads
#   3. discordant_reads
#   4. merge_evidence
#   5. annotate_candidates
# ============================================================

my $config;
my $bam_list;
my $outdir = "results";
my $force  = 0;
my $help   = 0;

GetOptions(
    "config=s"   => \$config,
    "bam-list=s" => \$bam_list,
    "outdir=s"   => \$outdir,
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

my %CONF = read_config($config);

my $perl     = get_conf(\%CONF, "PERL", "perl");
my $threads  = get_conf(\%CONF, "THREADS", 4);
my $keep_tmp = get_conf(\%CONF, "KEEP_TMP", 0);

my $exon_txt = resolve_path(
    get_conf(\%CONF, "REFSEQ_MANE_SELECT_EXON_TXT", ""),
    $Bin
);

my $gene_txt = resolve_path(
    get_conf(\%CONF, "REFSEQ_MANE_SELECT_GENE_TXT", ""),
    $Bin
);

die "[ERROR] Exon TXT file not found: $exon_txt\n" unless -s $exon_txt;
die "[ERROR] Gene TXT file not found: $gene_txt\n" unless -s $gene_txt;

check_exon_txt_format($exon_txt);
check_gene_txt_format($gene_txt);

$outdir = abs_path($outdir) || $outdir;
make_path($outdir) unless -d $outdir;
$outdir = abs_path($outdir);

my %SCRIPT = (
    depth      => "$Bin/bin/run_gene_mean_depth.pl",
    split      => "$Bin/bin/run_split_reads.pl",
    discordant => "$Bin/bin/run_discordant_reads.pl",
    merge      => "$Bin/bin/merge_evidence.pl",
    annotate   => "$Bin/bin/annotate_candidates.pl",
);

foreach my $key (sort keys %SCRIPT) {
    check_script($SCRIPT{$key});
}

my @samples = read_bam_list($bam_list);
die "[ERROR] No valid BAM found in BAM list: $bam_list\n" unless @samples;

my $shell_list = "$outdir/sample_shell.list";
my $qsub_list  = "$outdir/qsub_command.list";

open my $LIST, ">", $shell_list or die "[ERROR] Cannot write: $shell_list\n";
open my $QSUB, ">", $qsub_list  or die "[ERROR] Cannot write: $qsub_list\n";

print_summary(
    config     => $config,
    bam_list   => $bam_list,
    exon_txt   => $exon_txt,
    gene_txt   => $gene_txt,
    outdir     => $outdir,
    threads    => $threads,
    keep_tmp   => $keep_tmp,
    sample_num => scalar(@samples),
    shell_list => $shell_list,
    qsub_list  => $qsub_list,
);

foreach my $item (@samples) {
    my $sample = $item->{sample};
    my $bam    = $item->{bam};

    my %DIR = prepare_sample_dirs($outdir, $sample, $force);

    my %OUT = (
        depth      => "$DIR{depth}/$sample.depth_candidates.tsv",
        split      => "$DIR{split}/$sample.split_reads.tsv",
        discordant => "$DIR{discordant}/$sample.discordant_reads.tsv",
        merged     => "$DIR{candidate}/$sample.merged_candidates.tsv",
        annotated  => "$DIR{report}/$sample.annotated_candidates.tsv",
    );

    my $sample_shell = "$DIR{sample}/$sample.run.sh";

    write_sample_shell(
        shell      => $sample_shell,
        sample     => $sample,
        bam        => $bam,
        config     => $config,
        exon_txt   => $exon_txt,
        gene_txt   => $gene_txt,
        perl       => $perl,
        threads    => $threads,
        keep_tmp   => $keep_tmp,
        script_ref => \%SCRIPT,
        dir_ref    => \%DIR,
        out_ref    => \%OUT,
    );

    chmod 0755, $sample_shell
        or die "[ERROR] Cannot chmod sample shell: $sample_shell\n";

    print $LIST "$sample_shell\n";
    print $QSUB "qsub -cwd -l p=$threads,vf=4G -binding linear:$threads -N HCMExonDel_$sample $sample_shell\n";

    print "[INFO] Sample shell generated: $sample_shell\n";
}

close $LIST;
close $QSUB;

print "\n";
print "[INFO] All sample shells generated successfully\n";
print "[INFO] Shell list: $shell_list\n";
print "[INFO] Qsub list : $qsub_list\n";
print "\n";
print "Run locally:\n";
print "  while read sh; do bash \$sh; done < $shell_list\n";
print "\n";
print "Submit by qsub:\n";
print "  sh $qsub_list\n\n";

exit 0;

# ============================================================
# Main shell writer
# ============================================================

sub write_sample_shell {
    my %args = @_;

    my $shell      = $args{shell};
    my $sample     = $args{sample};
    my $bam        = $args{bam};
    my $config     = $args{config};
    my $exon_txt   = $args{exon_txt};
    my $gene_txt   = $args{gene_txt};
    my $perl       = $args{perl};
    my $threads    = $args{threads};
    my $keep_tmp   = $args{keep_tmp};
    my $script_ref = $args{script_ref};
    my $dir_ref    = $args{dir_ref};
    my $out_ref    = $args{out_ref};

    open my $SH, ">", $shell or die "[ERROR] Cannot write sample shell: $shell\n";

    write_header(
        fh       => $SH,
        sample   => $sample,
        config   => $config,
        bam      => $bam,
        exon_txt => $exon_txt,
        gene_txt => $gene_txt,
        outdir   => $dir_ref->{sample},
        threads  => $threads,
    );

    my $mkdir_cmd = join(" ",
        "mkdir -p",
        $dir_ref->{log},
        $dir_ref->{depth},
        $dir_ref->{split},
        $dir_ref->{discordant},
        $dir_ref->{candidate},
        $dir_ref->{report},
        $dir_ref->{tmp},
    );

    write_block($SH, "Start sample: $sample");
    print $SH "$mkdir_cmd\n\n";

    my $cmd;

    $cmd = "$perl $script_ref->{depth} "
         . "--config $config "
         . "--bam $bam "
         . "--sample $sample "
         . "--out $out_ref->{depth} "
         . "--threads $threads";
    write_cmd($SH, "Step 1: gene_mean_depth", $cmd, "$dir_ref->{log}/01.gene_mean_depth.log");

    $cmd = "$perl $script_ref->{split} "
         . "--config $config "
         . "--bam $bam "
         . "--sample $sample "
         . "--out $out_ref->{split} "
         . "--threads $threads";
    write_cmd($SH, "Step 2: split_reads", $cmd, "$dir_ref->{log}/02.split_reads.log");

    $cmd = "$perl $script_ref->{discordant} "
         . "--config $config "
         . "--bam $bam "
         . "--sample $sample "
         . "--out $out_ref->{discordant} "
         . "--threads $threads";
    write_cmd($SH, "Step 3: discordant_reads", $cmd, "$dir_ref->{log}/03.discordant_reads.log");

    $cmd = "$perl $script_ref->{merge} "
         . "--config $config "
         . "--sample $sample "
         . "--depth $out_ref->{depth} "
         . "--split $out_ref->{split} "
         . "--discordant $out_ref->{discordant} "
         . "--out $out_ref->{merged}";
    write_cmd($SH, "Step 4: merge_evidence", $cmd, "$dir_ref->{log}/04.merge_evidence.log");

    $cmd = "$perl $script_ref->{annotate} "
         . "--config $config "
         . "--input $out_ref->{merged} "
         . "--out $out_ref->{annotated}";
    write_cmd($SH, "Step 5: annotate_candidates", $cmd, "$dir_ref->{log}/05.annotate_candidates.log");

    if (!$keep_tmp) {
        print $SH "rm -rf $dir_ref->{tmp}\n\n";
    }

    print $SH "echo \"[INFO] Finished sample: $sample \$(date)\"\n";
    print $SH "echo \"[INFO] Annotated candidates: $out_ref->{annotated}\"\n";

    close $SH;
}

sub write_header {
    my %args = @_;

    my $FH      = $args{fh};
    my $sample  = $args{sample};
    my $config  = $args{config};
    my $bam     = $args{bam};
    my $exon    = $args{exon_txt};
    my $gene    = $args{gene_txt};
    my $outdir  = $args{outdir};
    my $threads = $args{threads};

    print $FH <<"HEADER";
#!/usr/bin/env bash
set -euo pipefail

############################################################
# Auto-generated by hcm_exondel.pl
# Sample   : $sample
# Config   : $config
# BAM      : $bam
# Exon TXT : $exon
# Gene TXT : $gene
# Outdir   : $outdir
# Threads  : $threads
############################################################

HEADER
}

sub write_block {
    my ($fh, $message) = @_;

    print $fh "echo \"[INFO] $message \$(date)\"\n";
}

sub write_cmd {
    my ($fh, $message, $cmd, $log) = @_;

    print $fh "echo \"[INFO] $message \$(date)\"\n";
    print $fh "$cmd > $log 2>&1\n\n";
}

# ============================================================
# Directory and input preparation
# ============================================================

sub prepare_sample_dirs {
    my ($outdir, $sample, $force) = @_;

    my %dir;

    $dir{sample}     = "$outdir/$sample";
    $dir{log}        = "$dir{sample}/00.log";
    $dir{depth}      = "$dir{sample}/01.depth";
    $dir{split}      = "$dir{sample}/02.split_reads";
    $dir{discordant} = "$dir{sample}/03.discordant_reads";
    $dir{candidate}  = "$dir{sample}/04.candidates";
    $dir{report}     = "$dir{sample}/05.report";
    $dir{tmp}        = "$dir{sample}/tmp";

    if (-d $dir{sample} && !$force) {
        die "[ERROR] Sample output directory already exists: $dir{sample}\n"
          . "        Please use --force to allow existing sample directories\n";
    }

    foreach my $key (keys %dir) {
        make_path($dir{$key}) unless -d $dir{$key};
    }

    return %dir;
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

# ============================================================
# Config and path functions
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

sub resolve_path {
    my ($path, $base_dir) = @_;

    die "[ERROR] Empty path provided\n"
        unless defined $path && $path ne "";

    if ($path =~ /^\//) {
        my $abs = abs_path($path);
        die "[ERROR] Path not found: $path\n" unless defined $abs;
        return $abs;
    }

    my $full = "$base_dir/$path";
    my $abs  = abs_path($full);

    die "[ERROR] Path not found: $full\n" unless defined $abs;

    return $abs;
}

# ============================================================
# Format checks
# ============================================================

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

# ============================================================
# Summary and usage
# ============================================================

sub print_summary {
    my %args = @_;

    print "\n";
    print "============================================================\n";
    print " HCMExonDel sample shell generation\n";
    print "============================================================\n";
    print "Config file : $args{config}\n";
    print "BAM list    : $args{bam_list}\n";
    print "Exon TXT    : $args{exon_txt}\n";
    print "Gene TXT    : $args{gene_txt}\n";
    print "Outdir      : $args{outdir}\n";
    print "Threads     : $args{threads}\n";
    print "Keep tmp    : $args{keep_tmp}\n";
    print "Samples     : $args{sample_num}\n";
    print "Shell list  : $args{shell_list}\n";
    print "Qsub list   : $args{qsub_list}\n";
    print "============================================================\n\n";
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
  --force        Overwrite existing sample directories.
  --help         Show this help message.

Output:
  <outdir>/sample_shell.list
  <outdir>/qsub_command.list
  <outdir>/<sample>/<sample>.run.sh

USAGE
}

