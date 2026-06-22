#!/usr/bin/env perl
# -*- coding: utf-8 -*-

use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);
use Cwd qw(abs_path);
use FindBin qw($Bin);
use File::Basename qw(dirname);

# ============================================================
# HCMExonDel main pipeline
#
# Function:
#   Generate one shell script for each BAM sample.
#
# Input:
#   1. Config file
#   2. BAM list, two TAB-delimited columns:
#      column 1: sample name
#      column 2: absolute BAM path
#
# Output:
#   OUTDIR/SAMPLE/SAMPLE.run.sh
#
# Workflow:
#   1. run_gene_mean_depth.pl
#   2. plot_depth_ratio.py
#   3. extract_sa_split_reads.pl
#   4. cluster_sa_split_reads.pl
#   5. run_discordant_reads.pl
#   6. merge_evidence.pl
#   7. extract_candidate_gene_bam.pl
#   8. annotate_candidates.pl
#
# Output structure:
#   OUTDIR/
#   └── SAMPLE/
#       ├── SAMPLE.run.sh
#       ├── 00.log/
#       ├── 01.depth/
#       ├── 02.split_reads/
#       ├── 03.discordant_reads/
#       ├── 04.candidates/
#       ├── 05.gene_bam/
#       ├── 06.report/
#       └── tmp/
#
# Notes:
#   1. No sample_shell.list is generated.
#   2. No qsub_command.list is generated.
#   3. THREADS is not read from config.
#   4. No --threads is passed to run_gene_mean_depth.pl.
#   5. Generated shell commands do not use single quotes around paths.
#   6. Candidate gene BAM extraction is always generated.
#   7. No flank is passed to extract_candidate_gene_bam.pl.
#      The extraction script uses its internal default flank = 0.
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
my $keep_tmp = normalize_bool(get_conf(\%CONF, "KEEP_TMP", 0));

my $project_root = detect_project_root($config);

my $exon_txt = resolve_config_path(
    get_conf_required(\%CONF, "REFSEQ_MANE_SELECT_EXON_TXT"),
    $project_root
);

my $gene_txt = resolve_config_path(
    get_conf_required(\%CONF, "REFSEQ_MANE_SELECT_GENE_TXT"),
    $project_root
);

die "[ERROR] REFSEQ_MANE_SELECT_EXON_TXT file not found: $exon_txt\n"
    unless -s $exon_txt;

die "[ERROR] REFSEQ_MANE_SELECT_GENE_TXT file not found: $gene_txt\n"
    unless -s $gene_txt;

check_exon_txt_format($exon_txt);
check_gene_txt_format($gene_txt);

$outdir = abs_path($outdir) || $outdir;
make_path($outdir) unless -d $outdir;
$outdir = abs_path($outdir);

my %SCRIPT = (
    depth         => "$Bin/bin/run_gene_mean_depth.pl",
    plot_depth    => "$Bin/bin/plot_depth_ratio.py",
    split_extract => "$Bin/bin/extract_sa_split_reads.pl",
    split_cluster => "$Bin/bin/cluster_sa_split_reads.pl",
    tlen_extract  => "$Bin/bin/extract_valid_tlen.pl",
    tlen_plot     => "$Bin/bin/plot_tlen_distribution.py",
    discordant    => "$Bin/bin/run_discordant_reads.pl",
    merge         => "$Bin/bin/merge_evidence.pl",
    extract_bam   => "$Bin/bin/extract_candidate_gene_bam.pl",
    annotate      => "$Bin/bin/annotate_candidates.pl",
);

foreach my $key (sort keys %SCRIPT) {
    check_script($SCRIPT{$key});
}

my @samples = read_bam_list($bam_list);
die "[ERROR] No valid sample found in BAM list: $bam_list\n" unless @samples;

print_summary(
    config       => $config,
    project_root => $project_root,
    bam_list     => $bam_list,
    exon_txt     => $exon_txt,
    gene_txt     => $gene_txt,
    outdir       => $outdir,
    keep_tmp     => $keep_tmp,
    sample_num   => scalar(@samples),
);

foreach my $item (@samples) {
    my $sample = $item->{sample};
    my $bam    = $item->{bam};

    my %DIR = prepare_sample_dirs($outdir, $sample, $force);

    my %OUT = (
        depth            => "$DIR{depth}/$sample.depth_candidates.tsv",
        all_window_ratio => "$DIR{depth}/$sample.depth_candidates.all_window_ratio.tsv",
        depth_plot_pdf   => "$DIR{depth}/$sample.window_del.per_gene.pdf",

        split_raw        => "$DIR{split}/$sample.split_reads.tsv",
        split_cluster    => "$DIR{split}/$sample.split_reads.clusters.tsv",

        valid_tlen       => "$DIR{discordant}/$sample.valid_tlen.tsv",
        tlen_plot_pdf    => "$DIR{discordant}/$sample.tlen_distribution.pdf",

        discordant       => "$DIR{discordant}/$sample.discordant_reads.tsv",

        merged           => "$DIR{candidate}/$sample.merged_candidates.tsv",
        annotated        => "$DIR{report}/$sample.annotated_candidates.tsv",
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
        keep_tmp   => $keep_tmp,
        script_ref => \%SCRIPT,
        dir_ref    => \%DIR,
        out_ref    => \%OUT,
    );

    chmod 0755, $sample_shell
        or die "[ERROR] Cannot chmod sample shell: $sample_shell\n";

    print "[INFO] Sample shell generated: $sample_shell\n";
}

print "\n";
print "[INFO] All sample shells generated successfully\n";
#print "[INFO] No sample_shell.list generated\n";
#print "[INFO] No qsub_command.list generated\n";
print "\n";
print "Run example:\n";
print "  bash $outdir/SAMPLE/SAMPLE.run.sh\n\n";

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
    my $keep_tmp   = $args{keep_tmp};
    my $script_ref = $args{script_ref};
    my $dir_ref    = $args{dir_ref};
    my $out_ref    = $args{out_ref};

    open my $SH, ">", $shell
        or die "[ERROR] Cannot write sample shell: $shell\n";

    write_header(
        fh       => $SH,
        sample   => $sample,
        config   => $config,
        bam      => $bam,
        exon_txt => $exon_txt,
        gene_txt => $gene_txt,
        outdir   => $dir_ref->{sample},
    );

    my $mkdir_cmd = join(
        " ",
        "mkdir -p",
        shell_quote($dir_ref->{log}),
        shell_quote($dir_ref->{depth}),
        shell_quote($dir_ref->{split}),
        shell_quote($dir_ref->{discordant}),
        shell_quote($dir_ref->{candidate}),
        shell_quote($dir_ref->{gene_bam}),
        shell_quote($dir_ref->{report}),
        shell_quote($dir_ref->{tmp}),
    );

    write_block($SH, "Start sample: $sample");
    print $SH "$mkdir_cmd\n\n";

    my $cmd;

    # ------------------------------------------------------------
    # Step 1. Depth candidates
    # ------------------------------------------------------------
    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{depth}),
        "--config", shell_quote($config),
        "--bam",    shell_quote($bam),
        "--sample", shell_quote($sample),
        "--out",    shell_quote($out_ref->{depth}),
    );

    write_cmd(
        $SH,
        "Step 1: gene_mean_depth",
        $cmd,
        "$dir_ref->{log}/01.gene_mean_depth.log"
    );

    # ------------------------------------------------------------
    # Step 1b. Plot depth ratio
    # ------------------------------------------------------------
    $cmd = join(
        " ",
        "python3",
        shell_quote($script_ref->{plot_depth}),
        "--conf",   shell_quote($config),
        "--input",  shell_quote($out_ref->{all_window_ratio}),
        "--output", shell_quote($out_ref->{depth_plot_pdf}),
    );

    write_cmd(
        $SH,
        "Step 1b: plot_depth_ratio",
        $cmd,
        "$dir_ref->{log}/01b.plot_depth_ratio.log"
    );

    # ------------------------------------------------------------
    # Step 2. Extract SA split reads
    # ------------------------------------------------------------
    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{split_extract}),
        "--conf",   shell_quote($config),
        "--bam",    shell_quote($bam),
        "--sample", shell_quote($sample),
        "--out",    shell_quote($out_ref->{split_raw}),
    );

    write_cmd(
        $SH,
        "Step 2: extract_sa_split_reads",
        $cmd,
        "$dir_ref->{log}/02.extract_sa_split_reads.log"
    );

    # ------------------------------------------------------------
    # Step 3. Cluster SA split reads
    # ------------------------------------------------------------
    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{split_cluster}),
        "--conf",    shell_quote($config),
        "--input",   shell_quote($out_ref->{split_raw}),
        "--outfile", shell_quote($out_ref->{split_cluster}),
    );

    write_cmd(
        $SH,
        "Step 3: cluster_sa_split_reads",
        $cmd,
        "$dir_ref->{log}/03.cluster_sa_split_reads.log"
    );

    # ------------------------------------------------------------
    # Step 3b. Extract valid TLEN
    # ------------------------------------------------------------
    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{tlen_extract}),
        shell_quote($bam),
        shell_quote($out_ref->{valid_tlen}),
    );

    write_cmd(
        $SH,
        "Step 3b: extract_valid_tlen",
        $cmd,
        "$dir_ref->{log}/03b.extract_valid_tlen.log"
    );

    # ------------------------------------------------------------
    # Step 3c. Plot TLEN distribution
    # ------------------------------------------------------------
    $cmd = join(
        " ",
        "python",
        shell_quote($script_ref->{tlen_plot}),
        shell_quote($out_ref->{valid_tlen}),
        shell_quote($out_ref->{tlen_plot_pdf}),
    );

    write_cmd(
        $SH,
        "Step 3c: plot_tlen_distribution",
        $cmd,
        "$dir_ref->{log}/03c.plot_tlen_distribution.log"
    );

    # ------------------------------------------------------------
    # Step 4. Discordant reads
    # ------------------------------------------------------------
    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{discordant}),
        "--config", shell_quote($config),
        "--bam",    shell_quote($bam),
        "--sample", shell_quote($sample),
        "--out",    shell_quote($out_ref->{discordant}),
    );

    write_cmd(
        $SH,
        "Step 4: discordant_reads",
        $cmd,
        "$dir_ref->{log}/04.discordant_reads.log"
    );

    # ------------------------------------------------------------
    # Step 5. Merge evidence
    # ------------------------------------------------------------
    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{merge}),
        "--config",     shell_quote($config),
        "--sample",     shell_quote($sample),
        "--depth",      shell_quote($out_ref->{depth}),
        "--split",      shell_quote($out_ref->{split_cluster}),
        "--discordant", shell_quote($out_ref->{discordant}),
        "--out",        shell_quote($out_ref->{merged}),
    );

    write_cmd(
        $SH,
        "Step 5: merge_evidence",
        $cmd,
        "$dir_ref->{log}/05.merge_evidence.log"
    );

    # ------------------------------------------------------------
    # Step 6. Extract candidate gene BAM
    # ------------------------------------------------------------
    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{extract_bam}),
        "--config",    shell_quote($config),
        "--candidate", shell_quote($out_ref->{merged}),
        "--bam",       shell_quote($bam),
        "--outdir",    shell_quote($dir_ref->{gene_bam}),
    );

    write_cmd(
        $SH,
        "Step 6: extract_candidate_gene_bam",
        $cmd,
        "$dir_ref->{log}/06.extract_candidate_gene_bam.log"
    );

    # ------------------------------------------------------------
    # Step 7. Annotate candidates
    # ------------------------------------------------------------
    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{annotate}),
        "--config", shell_quote($config),
        "--input",  shell_quote($out_ref->{merged}),
        "--out",    shell_quote($out_ref->{annotated}),
    );

    write_cmd(
        $SH,
        "Step 7: annotate_candidates",
        $cmd,
        "$dir_ref->{log}/07.annotate_candidates.log"
    );

    if (!$keep_tmp) {
        print $SH "rm -rf " . shell_quote($dir_ref->{tmp}) . "\n\n";
    }

    print $SH "echo \"[INFO] Finished sample: $sample \$(date)\"\n";
    print $SH "echo \"[INFO] Depth candidates       : $out_ref->{depth}\"\n";
    print $SH "echo \"[INFO] Raw split reads        : $out_ref->{split_raw}\"\n";
    print $SH "echo \"[INFO] Split clusters        : $out_ref->{split_cluster}\"\n";
    print $SH "echo \"[INFO] Discordant candidates : $out_ref->{discordant}\"\n";
    print $SH "echo \"[INFO] Merged candidates     : $out_ref->{merged}\"\n";
    print $SH "echo \"[INFO] Candidate gene BAM dir: $dir_ref->{gene_bam}\"\n";
    print $SH "echo \"[INFO] Annotated candidates  : $out_ref->{annotated}\"\n";

    close $SH;
}

sub write_header {
    my %args = @_;

    my $FH     = $args{fh};
    my $sample = $args{sample};
    my $config = $args{config};
    my $bam    = $args{bam};
    my $exon   = $args{exon_txt};
    my $gene   = $args{gene_txt};
    my $outdir = $args{outdir};

    print $FH <<"HEADER";
#!/usr/bin/env bash
set -euo pipefail

############################################################
# Auto-generated by HCMExonDel.pl
# Sample    : $sample
# Config    : $config
# BAM       : $bam
# Exon TXT  : $exon
# Gene TXT  : $gene
# Outdir    : $outdir
#
# Workflow:
#   1. gene_mean_depth
#   2. plot_depth_ratio
#   3. extract_sa_split_reads
#   4. cluster_sa_split_reads
#   5. discordant_reads
#   6. merge_evidence
#   7. extract_candidate_gene_bam
#   8. annotate_candidates
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
    print $fh "$cmd > " . shell_quote($log) . " 2>&1\n\n";
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
    $dir{gene_bam}   = "$dir{sample}/05.gene_bam";
    $dir{report}     = "$dir{sample}/06.report";
    $dir{tmp}        = "$dir{sample}/tmp";

    if (-d $dir{sample} && !$force) {
        die "[ERROR] Sample output directory already exists: $dir{sample}\n"
          . "        Please use --force to allow writing into existing sample directory\n";
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
    my %seen_bam;

    open my $FH, "<", $file
        or die "[ERROR] Cannot open BAM list: $file\n";

    my $line_no = 0;

    while (my $line = <$FH>) {
        chomp $line;
        $line =~ s/\r$//;
        $line_no++;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        my @fields = split /\t/, $line, -1;

        die "[ERROR] Invalid BAM list format at line $line_no:\n"
          . "        $line\n"
          . "        Expected two TAB-delimited columns:\n"
          . "        SampleName<TAB>/absolute/path/to/sample.bam\n"
          unless @fields == 2;

        my ($sample, $bam) = @fields;

        $sample =~ s/^\s+|\s+$//g;
        $bam    =~ s/^\s+|\s+$//g;

        die "[ERROR] Empty sample name at line $line_no\n"
            unless defined $sample && $sample ne "";

        die "[ERROR] Empty BAM path at line $line_no\n"
            unless defined $bam && $bam ne "";

        die "[ERROR] Invalid sample name at line $line_no: $sample\n"
          . "        Sample name can only contain letters, numbers, dot, underscore and hyphen.\n"
          unless $sample =~ /^[A-Za-z0-9_.-]+$/;

        die "[ERROR] BAM path must be an absolute path at line $line_no:\n"
          . "        $bam\n"
          unless $bam =~ m{^/};

        die "[ERROR] BAM file does not end with .bam at line $line_no:\n"
          . "        $bam\n"
          unless $bam =~ /\.bam$/i;

        my $abs_bam = abs_path($bam);

        die "[ERROR] BAM file not found at line $line_no:\n"
          . "        $bam\n"
          unless $abs_bam && -s $abs_bam;

        check_bam_index($abs_bam);

        die "[ERROR] Duplicated sample name detected at line $line_no: $sample\n"
            if exists $seen_sample{$sample};

        die "[ERROR] Duplicated BAM path detected at line $line_no:\n"
          . "        $abs_bam\n"
            if exists $seen_bam{$abs_bam};

        $seen_sample{$sample} = 1;
        $seen_bam{$abs_bam}   = 1;

        push @list, {
            sample => $sample,
            bam    => $abs_bam,
        };
    }

    close $FH;

    return @list;
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

    open my $FH, "<", $file
        or die "[ERROR] Cannot open config file: $file\n";

    while (my $line = <$FH>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        $line =~ s/\s+#.*$//;

        if ($line =~ /^\s*([^=\s]+)\s*=\s*(.*?)\s*$/) {
            my $key = $1;
            my $val = $2;

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

    if (
        exists $conf_ref->{$key}
        && defined $conf_ref->{$key}
        && $conf_ref->{$key} ne ""
    ) {
        return $conf_ref->{$key};
    }

    return $default;
}

sub get_conf_required {
    my ($conf_ref, $key) = @_;

    die "[ERROR] Required config parameter missing: $key\n"
        unless exists $conf_ref->{$key}
        && defined $conf_ref->{$key}
        && $conf_ref->{$key} ne "";

    return $conf_ref->{$key};
}

sub normalize_bool {
    my ($v) = @_;

    return 0 unless defined $v;

    $v =~ s/^\s+|\s+$//g;

    return 1 if $v =~ /^(1|yes|true|on)$/i;
    return 0 if $v =~ /^(0|no|false|off)$/i;

    return $v ? 1 : 0;
}

sub detect_project_root {
    my ($config_path) = @_;

    my $config_dir = dirname($config_path);

    if ($config_dir =~ m{/conf$}) {
        return abs_path(dirname($config_dir));
    }

    return abs_path($config_dir);
}

sub resolve_config_path {
    my ($path, $project_root) = @_;

    die "[ERROR] Empty path provided\n"
        unless defined $path && $path ne "";

    if ($path =~ m{^/}) {
        my $abs = abs_path($path);
        die "[ERROR] Path not found: $path\n" unless defined $abs;
        return $abs;
    }

    my $full = "$project_root/$path";
    my $abs  = abs_path($full);

    die "[ERROR] Path not found: $full\n" unless defined $abs;

    return $abs;
}

sub shell_quote {
    my ($str) = @_;

    die "[ERROR] Undefined shell argument\n" unless defined $str;

    $str =~ s/^\s+|\s+$//g;

    die "[ERROR] Empty shell argument\n" if $str eq "";

    # Do not wrap arguments with single quotes.
    # Escape common shell metacharacters only when needed.
    $str =~ s/([ \t\n\r\\\"\`\$\&\|\;\<\>\(\)\{\}\[\]\*\?\!\#])/\\$1/g;

    return $str;
}

# ============================================================
# Checks
# ============================================================

sub check_script {
    my ($script) = @_;

    die "[ERROR] Required script not found: $script\n"
        unless -s $script;

    return 1;
}

sub check_exon_txt_format {
    my ($file) = @_;

    open my $FH, "<", $file
        or die "[ERROR] Cannot open exon TXT: $file\n";

    my $header = <$FH>;
    close $FH;

    die "[ERROR] Empty exon TXT file: $file\n" unless defined $header;

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

    open my $FH, "<", $file
        or die "[ERROR] Cannot open gene TXT: $file\n";

    my $header = <$FH>;
    close $FH;

    die "[ERROR] Empty gene TXT file: $file\n" unless defined $header;

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
    print "Config file              : $args{config}\n";
    print "Project root             : $args{project_root}\n";
    print "BAM list                 : $args{bam_list}\n";
    print "BAM list format          : SampleName<TAB>/absolute/path/to/sample.bam\n";
    print "Exon TXT                 : $args{exon_txt}\n";
    print "Gene TXT                 : $args{gene_txt}\n";
    print "Outdir                   : $args{outdir}\n";
    print "Keep tmp                 : $args{keep_tmp}\n";
    print "Samples                  : $args{sample_num}\n";
    print "Workflow                 : depth -> plot -> SA split -> SA cluster -> discordant -> merge -> gene BAM -> annotate\n";
    print "Generated output         : OUTDIR/SAMPLE/SAMPLE.run.sh\n";
    #print "sample_shell.list        : not generated\n";
    #print "qsub_command.list        : not generated\n";
    #print "Split input for merge    : *.split_reads.clusters.tsv\n";
    #print "Candidate gene BAM       : enabled\n";
    #print "Candidate gene BAM flank : default 0, strict gene coordinates\n";
    #print "Shell quoting style      : no single quotes, escape metacharacters only\n";
    print "============================================================\n\n";
}

sub usage {
    return <<"USAGE";
Usage:
  perl HCMExonDel.pl \\
    --config conf/hcm_exondel.example.conf \\
    --bam-list input_bam.list \\
    --outdir results \\
    --force

Required arguments:
  --config      Configuration file.
  --bam-list    BAM list file. Two TAB-delimited columns:
                column 1: sample name
                column 2: absolute BAM path

Optional arguments:
  --outdir      Output directory. [default: results]
  --force       Allow writing into existing sample directories.
  --help        Show this help message.

Required config keys checked by HCMExonDel.pl:
  REFSEQ_MANE_SELECT_EXON_TXT
  REFSEQ_MANE_SELECT_GENE_TXT

Config keys used directly by HCMExonDel.pl:
  PERL
  KEEP_TMP

Other config keys are checked by each bin script itself, for example:
  SAMTOOLS

  HCM_CORE_GENE_LIST
  ANALYZE_CORE_GENES_ONLY
  TARGET_REGION_FLANK

  MIN_MAPQ
  EXCLUDE_DUPLICATES

  WINDOW_SIZE
  WINDOW_STEP
  MIN_WINDOW_SIZE
  MIN_GENE_MEAN_DEPTH
  DEL_DEPTH_RATIO_CUTOFF
  MIN_CONSECUTIVE_DEL_WINDOWS
  KEEP_GENE_DEPTH_FILE

  SA_SPLIT_CLUSTER_WINDOW
  SA_SPLIT_MIN_SUPPORT_READS

  MIN_DISCORDANT_INSERT_SIZE
  MIN_DISCORDANT_READS
  DISCORDANT_CLUSTER_DISTANCE
  FILTER_DELETION_ORIENTATION
  REF_FASTA_INDEX

  MIN_EXON_OVERLAP_FRACTION

  VERBOSE

BAM list format:
  SampleName<TAB>/absolute/path/to/sample.bam

Example:
  25B09089386<TAB>/ehpcdata/fulongfei/project/XJ_HCM_WGS_FHOD3/JX_2/25B09089386.final.merge.bam

Workflow:
  1. run_gene_mean_depth.pl
  2. plot_depth_ratio.py
  3. extract_sa_split_reads.pl
  4. cluster_sa_split_reads.pl
  5. run_discordant_reads.pl
  6. merge_evidence.pl
  7. extract_candidate_gene_bam.pl
  8. annotate_candidates.pl

Output:
  OUTDIR/SAMPLE/SAMPLE.run.sh

Output structure:
  OUTDIR/
  └── SAMPLE/
      ├── SAMPLE.run.sh
      ├── 00.log/
      ├── 01.depth/
      │   ├── SAMPLE.depth_candidates.tsv
      │   ├── SAMPLE.depth_candidates.all_window_ratio.tsv
      │   ├── SAMPLE.depth_candidates.del_windows.tsv
      │   ├── SAMPLE.depth_candidates.gene_mean_depth.tsv
      │   ├── SAMPLE.depth_candidates.window_depth.tsv
      │   ├── SAMPLE.depth_candidates.depth_ratio.pdf
      │   └── SAMPLE.depth_candidates.gene_depth_files/
      ├── 02.split_reads/
      │   ├── SAMPLE.split_reads.tsv
      │   ├── SAMPLE.split_reads.clusters.tsv
      │   ├── SAMPLE.split_reads.failed_clusters.tsv
      │   └── SAMPLE.split_reads.supporting_reads.tsv
      ├── 03.discordant_reads/
      │   ├── SAMPLE.discordant_reads.tsv
      │   ├── SAMPLE.discordant_reads.raw_discordant_pairs.tsv
      │   ├── SAMPLE.discordant_reads.supporting_pairs.tsv
      │   └── SAMPLE.discordant_reads.discarded_clusters.tsv
      ├── 04.candidates/
      │   ├── SAMPLE.merged_candidates.tsv
      │   └── SAMPLE.merged_candidates.tsv.all.tsv
      ├── 05.gene_bam/
      │   ├── GENE.bam
      │   └── GENE.bam.bai
      ├── 06.report/
      │   └── SAMPLE.annotated_candidates.tsv
      └── tmp/

Important:
  No sample_shell.list is generated.
  No qsub_command.list is generated.

  THREADS is not used by HCMExonDel.pl.
  No --threads is passed to run_gene_mean_depth.pl.

  merge_evidence.pl uses:
    --split 02.split_reads/SAMPLE.split_reads.clusters.tsv

  extract_candidate_gene_bam.pl is always generated.
  No --flank is passed, so it uses its own default flank = 0.
  BAM extraction is strictly based on gene coordinates.

  Generated shell commands do not wrap paths with single quotes.
  Common shell metacharacters are escaped with backslashes when needed.

USAGE
}

