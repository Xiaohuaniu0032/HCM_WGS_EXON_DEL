#!/usr/bin/env perl
# -*- coding: utf-8 -*-

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use FindBin qw($Bin);

# ============================================================
# HCMExonDel.pl
#
# Main pipeline generator for HCM exon-level deletion detection.
#
# Function:
#   Generate one executable run.sh for each BAM sample.
#
# Input:
#   1. Config file
#   2. BAM list, two TAB-delimited columns:
#        SampleID    /absolute/path/to/sample.bam
#
# Output:
#   OUTDIR/SAMPLE/SAMPLE.run.sh
#
# Workflow:
#   Step 00. extract_target_bam
#   Step 01. gene_mean_depth
#   Step 02. cusum_depth_evidence
#   Step 03. base_depth_ratio
#   Step 04. plot_depth_ratio
#   Step 05. extract_sa_split_reads
#   Step 06. cluster_sa_split_reads
#   Step 07. extract_valid_tlen
#   Step 08. plot_tlen_distribution
#   Step 09. discordant_reads
#   Step 10. merge_evidence
#   Step 11. annotate_candidates
#
# Design notes:
#   1. HCMExonDel only analyzes HCM core genes.
#   2. The original WGS BAM is used only by extract_target_bam.pl.
#   3. All downstream analyses use SAMPLE.target.bam.
#   4. merge_evidence.pl is split-read-centered:
#        split cluster is used as the candidate event backbone;
#        depth and discordant signals are used as supporting evidence.
#   5. merge_evidence.pl does not read config directly.
#      HCMExonDel.pl reads EVIDENCE_OVERLAP_FRACTION from config
#      and passes it to merge_evidence.pl as --min-overlap.
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
    "help|h"     => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless defined $config && defined $bam_list;

$config = abs_path($config);
die "[ERROR] Config file not found: $config\n"
    unless defined $config && -s $config;

$bam_list = abs_path($bam_list);
die "[ERROR] BAM list file not found: $bam_list\n"
    unless defined $bam_list && -s $bam_list;

my %CONF = read_config($config);

# ------------------------------------------------------------
# Runtime executables
# ------------------------------------------------------------
my $perl   = get_conf_required(\%CONF, "PERL");
my $python = get_conf_required(\%CONF, "PYTHON");

check_executable($perl,   "PERL");
check_executable($python, "PYTHON");

# ------------------------------------------------------------
# Workflow switches
# ------------------------------------------------------------
my $keep_tmp = normalize_bool(
    get_conf_required(\%CONF, "KEEP_TMP"),
    "KEEP_TMP"
);

my $keep_gene_depth_file = normalize_bool(
    get_conf_required(\%CONF, "KEEP_GENE_DEPTH_FILE"),
    "KEEP_GENE_DEPTH_FILE"
);

my $output_base_depth_ratio = normalize_bool(
    get_conf_required(\%CONF, "OUTPUT_BASE_DEPTH_RATIO"),
    "OUTPUT_BASE_DEPTH_RATIO"
);

my $output_base_depth_ratio_all = normalize_bool(
    get_conf_required(\%CONF, "OUTPUT_BASE_DEPTH_RATIO_ALL"),
    "OUTPUT_BASE_DEPTH_RATIO_ALL"
);

my $run_tlen_qc = normalize_bool(
    get_conf_required(\%CONF, "RUN_TLEN_QC"),
    "RUN_TLEN_QC"
);

die "[ERROR] KEEP_GENE_DEPTH_FILE must be 1 because Step 02 CUSUM depth evidence depends on per-gene depth files.\n"
    unless $keep_gene_depth_file;

die "[ERROR] OUTPUT_BASE_DEPTH_RATIO_ALL=1 requires OUTPUT_BASE_DEPTH_RATIO=1.\n"
    if $output_base_depth_ratio_all && !$output_base_depth_ratio;

# ------------------------------------------------------------
# Deprecated parameters
# ------------------------------------------------------------
if (exists $CONF{ANALYZE_CORE_GENES_ONLY}) {
    die "[ERROR] ANALYZE_CORE_GENES_ONLY has been removed.\n"
      . "        HCMExonDel now always analyzes HCM core genes from HCM_CORE_GENE_LIST.\n"
      . "        Please delete ANALYZE_CORE_GENES_ONLY from the config file.\n";
}

if (exists $CONF{CANDIDATE_BAM_FLANK}) {
    die "[ERROR] CANDIDATE_BAM_FLANK has been removed because candidate gene BAM extraction is no longer part of the workflow.\n"
      . "        Please delete CANDIDATE_BAM_FLANK from the config file.\n";
}

# ------------------------------------------------------------
# Algorithm parameters checked by main pipeline
# ------------------------------------------------------------
my $target_region_flank = get_conf_required(\%CONF, "TARGET_REGION_FLANK");
validate_nonnegative_integer("TARGET_REGION_FLANK", $target_region_flank);

my $target_bam_threads = get_conf_required(\%CONF, "TARGET_BAM_THREADS");
validate_positive_integer("TARGET_BAM_THREADS", $target_bam_threads);

my $split_read_threads = get_conf_required(\%CONF, "SPLIT_READ_THREADS");
validate_positive_integer("SPLIT_READ_THREADS", $split_read_threads);

my $discordant_read_threads = get_conf_required(\%CONF, "DISCORDANT_READ_THREADS");
validate_positive_integer("DISCORDANT_READ_THREADS", $discordant_read_threads);

my $evidence_overlap_fraction = get_conf_required(\%CONF, "EVIDENCE_OVERLAP_FRACTION");
validate_fraction("EVIDENCE_OVERLAP_FRACTION", $evidence_overlap_fraction);

my %CUSUM = (
    baseline      => get_conf_required(\%CONF, "CUSUM_BASELINE"),
    bin_size      => get_conf_required(\%CONF, "CUSUM_BIN_SIZE"),
    k             => get_conf_required(\%CONF, "CUSUM_K"),
    h             => get_conf_required(\%CONF, "CUSUM_H"),
    del_ratio     => get_conf_required(\%CONF, "CUSUM_DEL_RATIO"),
    min_bins      => get_conf_required(\%CONF, "CUSUM_MIN_BINS"),
    min_len       => get_conf_required(\%CONF, "CUSUM_MIN_LEN"),
    edge_ratio    => get_conf_required(\%CONF, "CUSUM_EDGE_RATIO"),
    recover_ratio => get_conf_required(\%CONF, "CUSUM_RECOVER_RATIO"),
    recover_bins  => get_conf_required(\%CONF, "CUSUM_RECOVER_BINS"),
);

validate_cusum_params(\%CUSUM);

my %TLEN;

if ($run_tlen_qc) {
    %TLEN = (
        max_pairs   => get_conf_required(\%CONF, "TLEN_MAX_PAIRS"),
        threads     => get_conf_required(\%CONF, "TLEN_THREADS"),
        random_seed => get_conf_required(\%CONF, "TLEN_RANDOM_SEED"),
    );

    validate_positive_integer("TLEN_MAX_PAIRS",  $TLEN{max_pairs});
    validate_positive_integer("TLEN_THREADS",    $TLEN{threads});
    validate_integer("TLEN_RANDOM_SEED",         $TLEN{random_seed});
}

# ------------------------------------------------------------
# Annotation files
# ------------------------------------------------------------
my $project_root = detect_project_root($config);

my $exon_txt = resolve_config_path(
    get_conf_required(\%CONF, "REFSEQ_MANE_SELECT_EXON_TXT"),
    $project_root
);

my $gene_txt = resolve_config_path(
    get_conf_required(\%CONF, "REFSEQ_MANE_SELECT_GENE_TXT"),
    $project_root
);

my $core_gene_list = resolve_config_path(
    get_conf_required(\%CONF, "HCM_CORE_GENE_LIST"),
    $project_root
);

die "[ERROR] REFSEQ_MANE_SELECT_EXON_TXT file not found: $exon_txt\n"
    unless -s $exon_txt;

die "[ERROR] REFSEQ_MANE_SELECT_GENE_TXT file not found: $gene_txt\n"
    unless -s $gene_txt;

die "[ERROR] HCM_CORE_GENE_LIST file not found: $core_gene_list\n"
    unless -s $core_gene_list;

check_exon_txt_format($exon_txt);
check_gene_txt_format($gene_txt);
check_core_gene_list_format($core_gene_list);

# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------
$outdir = abs_path($outdir) || $outdir;

make_path($outdir) unless -d $outdir;

$outdir = abs_path($outdir);

# ------------------------------------------------------------
# Pipeline scripts
# ------------------------------------------------------------
my %SCRIPT = (
    target_bam        => "$Bin/bin/extract_target_bam.pl",
    depth             => "$Bin/bin/run_gene_mean_depth.pl",
    sample_cusum      => "$Bin/bin/run_sample_cusum_depth.pl",
    cusum_depth       => "$Bin/bin/cusum_depth_del.pl",
    sample_base_ratio => "$Bin/bin/run_sample_base_depth_ratio.pl",
    calc_ratio        => "$Bin/bin/calc_depth_ratio.pl",
    plot_depth        => "$Bin/bin/plot_depth_ratio.py",
    split_extract     => "$Bin/bin/extract_sa_split_reads.pl",
    split_cluster     => "$Bin/bin/cluster_sa_split_reads.pl",
    tlen_extract      => "$Bin/bin/extract_valid_tlen.pl",
    tlen_plot         => "$Bin/bin/plot_tlen_distribution.py",
    discordant        => "$Bin/bin/run_discordant_reads.pl",
    merge             => "$Bin/bin/merge_evidence.pl",
    annotate          => "$Bin/bin/annotate_candidates.pl",
);

my @required_scripts = qw/
    target_bam
    depth
    sample_cusum
    cusum_depth
    plot_depth
    split_extract
    split_cluster
    discordant
    merge
    annotate
/;

push @required_scripts, qw/sample_base_ratio calc_ratio/
    if $output_base_depth_ratio;

push @required_scripts, qw/tlen_extract tlen_plot/
    if $run_tlen_qc;

for my $key (@required_scripts) {
    check_script($SCRIPT{$key}, $key);
}

# ------------------------------------------------------------
# BAM list
# ------------------------------------------------------------
my @samples = read_bam_list($bam_list);

die "[ERROR] No valid sample found in BAM list: $bam_list\n"
    unless @samples;

print_summary(
    config                      => $config,
    project_root                => $project_root,
    bam_list                    => $bam_list,
    exon_txt                    => $exon_txt,
    gene_txt                    => $gene_txt,
    core_gene_list              => $core_gene_list,
    outdir                      => $outdir,
    perl                        => $perl,
    python                      => $python,
    keep_tmp                    => $keep_tmp,
    keep_gene_depth_file        => $keep_gene_depth_file,
    output_base_depth_ratio     => $output_base_depth_ratio,
    output_base_depth_ratio_all => $output_base_depth_ratio_all,
    run_tlen_qc                 => $run_tlen_qc,
    target_region_flank         => $target_region_flank,
    target_bam_threads          => $target_bam_threads,
    split_read_threads          => $split_read_threads,
    discordant_read_threads     => $discordant_read_threads,
    evidence_overlap_fraction   => $evidence_overlap_fraction,
    cusum_ref                   => \%CUSUM,
    tlen_ref                    => \%TLEN,
    sample_num                  => scalar(@samples),
);

for my $item (@samples) {
    my $sample = $item->{sample};
    my $bam    = $item->{bam};

    my %DIR = prepare_sample_dirs($outdir, $sample, $force);
    my %OUT = build_sample_outputs($sample, \%DIR);

    my $sample_shell = "$DIR{sample}/$sample.run.sh";

    write_sample_shell(
        shell                       => $sample_shell,
        sample                      => $sample,
        bam                         => $bam,
        config                      => $config,
        exon_txt                    => $exon_txt,
        gene_txt                    => $gene_txt,
        perl                        => $perl,
        python                      => $python,
        keep_tmp                    => $keep_tmp,
        output_base_depth_ratio     => $output_base_depth_ratio,
        output_base_depth_ratio_all => $output_base_depth_ratio_all,
        run_tlen_qc                 => $run_tlen_qc,
        evidence_overlap_fraction   => $evidence_overlap_fraction,
        script_ref                  => \%SCRIPT,
        dir_ref                     => \%DIR,
        out_ref                     => \%OUT,
    );

    chmod 0755, $sample_shell
        or die "[ERROR] Cannot chmod sample shell: $sample_shell\n";

    print "[INFO] Sample shell generated: $sample_shell\n";
}

print "\n";
print "[INFO] All sample shells generated successfully\n\n";
print "Run example:\n";
print "  bash $outdir/SAMPLE/SAMPLE.run.sh\n\n";

exit 0;

# ============================================================
# Main shell writer
# ============================================================

sub write_sample_shell {
    my %args = @_;

    my $shell                       = $args{shell};
    my $sample                      = $args{sample};
    my $bam                         = $args{bam};
    my $config                      = $args{config};
    my $exon_txt                    = $args{exon_txt};
    my $gene_txt                    = $args{gene_txt};
    my $perl                        = $args{perl};
    my $python                      = $args{python};
    my $keep_tmp                    = $args{keep_tmp};
    my $output_base_depth_ratio     = $args{output_base_depth_ratio};
    my $output_base_depth_ratio_all = $args{output_base_depth_ratio_all};
    my $run_tlen_qc                 = $args{run_tlen_qc};
    my $evidence_overlap_fraction   = $args{evidence_overlap_fraction};
    my $script_ref                  = $args{script_ref};
    my $dir_ref                     = $args{dir_ref};
    my $out_ref                     = $args{out_ref};

    open my $SH, ">", $shell
        or die "[ERROR] Cannot write sample shell: $shell\n";

    write_header(
        fh           => $SH,
        sample       => $sample,
        config       => $config,
        original_bam => $bam,
        target_bam   => $out_ref->{target_bam},
        exon_txt     => $exon_txt,
        gene_txt     => $gene_txt,
        outdir       => $dir_ref->{sample},
    );

    my $mkdir_cmd = join(
        " ",
        "mkdir -p",
        shell_quote($dir_ref->{log}),
        shell_quote($dir_ref->{target_bam}),
        shell_quote($dir_ref->{depth}),
        shell_quote($out_ref->{cusum_dir}),
        shell_quote($out_ref->{base_depth_ratio_dir}),
        shell_quote($dir_ref->{split}),
        shell_quote($dir_ref->{discordant}),
        shell_quote($dir_ref->{candidate}),
        shell_quote($dir_ref->{report}),
        shell_quote($dir_ref->{tmp}),
    );

    write_block($SH, "Start sample: $sample");
    print $SH "$mkdir_cmd\n\n";

    my $cmd;

    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{target_bam}),
        "--config", shell_quote($config),
        "--sample", shell_quote($sample),
        "--bam",    shell_quote($bam),
        "--outdir", shell_quote($dir_ref->{target_bam}),
    );
    write_cmd(
        $SH,
        "Step 00: extract_target_bam",
        $cmd,
        "$dir_ref->{log}/00.extract_target_bam.log"
    );

    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{depth}),
        "--config", shell_quote($config),
        "--bam",    shell_quote($out_ref->{target_bam}),
        "--sample", shell_quote($sample),
        "--out",    shell_quote($out_ref->{depth}),
    );
    write_cmd(
        $SH,
        "Step 01: gene_mean_depth",
        $cmd,
        "$dir_ref->{log}/01.gene_mean_depth.log"
    );

    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{sample_cusum}),
        "--config",  shell_quote($config),
        "--in-dir",  shell_quote($out_ref->{gene_depth_dir}),
        "--out-dir", shell_quote($out_ref->{cusum_dir}),
        "--out",     shell_quote($out_ref->{cusum_all}),
        "--perl",    shell_quote($perl),
        "--script",  shell_quote($script_ref->{cusum_depth}),
    );
    write_cmd(
        $SH,
        "Step 02: cusum_depth_evidence",
        $cmd,
        "$dir_ref->{log}/02.cusum_depth_evidence.log"
    );

    if ($output_base_depth_ratio) {
        $cmd = join(
            " ",
            shell_quote($perl),
            shell_quote($script_ref->{sample_base_ratio}),
            "--in-dir",  shell_quote($out_ref->{gene_depth_dir}),
            "--out-dir", shell_quote($out_ref->{base_depth_ratio_dir}),
            "--perl",    shell_quote($perl),
            "--script",  shell_quote($script_ref->{calc_ratio}),
        );

        if ($output_base_depth_ratio_all) {
            $cmd .= " --ratio-all "   . shell_quote($out_ref->{base_depth_ratio_all});
            $cmd .= " --summary-all " . shell_quote($out_ref->{base_depth_summary_all});
        }

        write_cmd(
            $SH,
            "Step 03: base_depth_ratio",
            $cmd,
            "$dir_ref->{log}/03.base_depth_ratio.log"
        );
    }
    else {
        $cmd = "echo " . shell_quote("[INFO] OUTPUT_BASE_DEPTH_RATIO=0, skip Step 03 base_depth_ratio");

        write_cmd(
            $SH,
            "Step 03: base_depth_ratio skipped",
            $cmd,
            "$dir_ref->{log}/03.base_depth_ratio.log"
        );
    }

    $cmd = join(
        " ",
        shell_quote($python),
        shell_quote($script_ref->{plot_depth}),
        "--conf",   shell_quote($config),
        "--input",  shell_quote($out_ref->{all_window_ratio}),
        "--output", shell_quote($out_ref->{depth_plot_pdf}),
    );
    write_cmd(
        $SH,
        "Step 04: plot_depth_ratio",
        $cmd,
        "$dir_ref->{log}/04.plot_depth_ratio.log"
    );

    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{split_extract}),
        "--conf",   shell_quote($config),
        "--bam",    shell_quote($out_ref->{target_bam}),
        "--sample", shell_quote($sample),
        "--out",    shell_quote($out_ref->{split_raw}),
    );
    write_cmd(
        $SH,
        "Step 05: extract_sa_split_reads",
        $cmd,
        "$dir_ref->{log}/05.extract_sa_split_reads.log"
    );

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
        "Step 06: cluster_sa_split_reads",
        $cmd,
        "$dir_ref->{log}/06.cluster_sa_split_reads.log"
    );

    if ($run_tlen_qc) {
        $cmd = join(
            " ",
            shell_quote($perl),
            shell_quote($script_ref->{tlen_extract}),
            "--config", shell_quote($config),
            "--bam",    shell_quote($out_ref->{target_bam}),
            "--out",    shell_quote($out_ref->{valid_tlen}),
        );
        write_cmd(
            $SH,
            "Step 07: extract_valid_tlen",
            $cmd,
            "$dir_ref->{log}/07.extract_valid_tlen.log"
        );
    }
    else {
        $cmd = "echo " . shell_quote("[INFO] RUN_TLEN_QC=0, skip Step 07 extract_valid_tlen");

        write_cmd(
            $SH,
            "Step 07: extract_valid_tlen skipped",
            $cmd,
            "$dir_ref->{log}/07.extract_valid_tlen.log"
        );
    }

    if ($run_tlen_qc) {
        $cmd = join(
            " ",
            shell_quote($python),
            shell_quote($script_ref->{tlen_plot}),
            shell_quote($out_ref->{valid_tlen}),
            shell_quote($out_ref->{tlen_plot_pdf}),
        );
        write_cmd(
            $SH,
            "Step 08: plot_tlen_distribution",
            $cmd,
            "$dir_ref->{log}/08.plot_tlen_distribution.log"
        );
    }
    else {
        $cmd = "echo " . shell_quote("[INFO] RUN_TLEN_QC=0, skip Step 08 plot_tlen_distribution");

        write_cmd(
            $SH,
            "Step 08: plot_tlen_distribution skipped",
            $cmd,
            "$dir_ref->{log}/08.plot_tlen_distribution.log"
        );
    }

    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{discordant}),
        "--config", shell_quote($config),
        "--bam",    shell_quote($out_ref->{target_bam}),
        "--sample", shell_quote($sample),
        "--out",    shell_quote($out_ref->{discordant}),
    );
    write_cmd(
        $SH,
        "Step 09: discordant_reads",
        $cmd,
        "$dir_ref->{log}/09.discordant_reads.log"
    );

    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{merge}),
        "--sample",      shell_quote($sample),
        "--depth",       shell_quote($out_ref->{cusum_all}),
        "--split",       shell_quote($out_ref->{split_cluster}),
        "--discordant",  shell_quote($out_ref->{discordant}),
        "--min-overlap", shell_quote($evidence_overlap_fraction),
        "--out",         shell_quote($out_ref->{merged}),
    );
    write_cmd(
        $SH,
        "Step 10: merge_evidence",
        $cmd,
        "$dir_ref->{log}/10.merge_evidence.log"
    );

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
        "Step 11: annotate_candidates",
        $cmd,
        "$dir_ref->{log}/11.annotate_candidates.log"
    );

    if (!$keep_tmp) {
        print $SH "rm -rf " . shell_quote($dir_ref->{tmp}) . "\n\n";
    }

    print $SH "echo \"[INFO] Finished sample: $sample \$(date)\"\n";
    print $SH "echo \"[INFO] Original BAM : $bam\"\n";
    print $SH "echo \"[INFO] Target BAM : $out_ref->{target_bam}\"\n";
    print $SH "echo \"[INFO] Target BED : $out_ref->{target_bed}\"\n";
    print $SH "echo \"[INFO] Window depth candidates : $out_ref->{depth}\"\n";
    print $SH "echo \"[INFO] CUSUM depth evidence : $out_ref->{cusum_all}\"\n";
    print $SH "echo \"[INFO] Depth ratio plot : $out_ref->{depth_plot_pdf}\"\n";

    if ($output_base_depth_ratio) {
        print $SH "echo \"[INFO] Base depth ratio dir : $out_ref->{base_depth_ratio_dir}\"\n";

        if ($output_base_depth_ratio_all) {
            print $SH "echo \"[INFO] Base depth ratio all : $out_ref->{base_depth_ratio_all}\"\n";
            print $SH "echo \"[INFO] Base depth summary all : $out_ref->{base_depth_summary_all}\"\n";
        }
    }

    print $SH "echo \"[INFO] Raw split reads : $out_ref->{split_raw}\"\n";
    print $SH "echo \"[INFO] Split clusters : $out_ref->{split_cluster}\"\n";

    if ($run_tlen_qc) {
        print $SH "echo \"[INFO] Valid TLEN : $out_ref->{valid_tlen}\"\n";
        print $SH "echo \"[INFO] TLEN plot : $out_ref->{tlen_plot_pdf}\"\n";
    }

    print $SH "echo \"[INFO] Discordant candidates : $out_ref->{discordant}\"\n";
    print $SH "echo \"[INFO] Merged candidates : $out_ref->{merged}\"\n";
    print $SH "echo \"[INFO] Annotated candidates : $out_ref->{annotated}\"\n";

    close $SH;
}

sub write_header {
    my %args = @_;

    my $FH           = $args{fh};
    my $sample       = $args{sample};
    my $config       = $args{config};
    my $original_bam = $args{original_bam};
    my $target_bam   = $args{target_bam};
    my $exon_txt     = $args{exon_txt};
    my $gene_txt     = $args{gene_txt};
    my $outdir       = $args{outdir};

    print $FH <<"HEADER";
#!/usr/bin/env bash
set -euo pipefail

############################################################
# Auto-generated by HCMExonDel.pl
#
# Sample       : $sample
# Config       : $config
# Original BAM : $original_bam
# Target BAM   : $target_bam
# Exon TXT     : $exon_txt
# Gene TXT     : $gene_txt
# Outdir       : $outdir
#
# Workflow:
#   Step 00. extract_target_bam
#   Step 01. gene_mean_depth
#   Step 02. cusum_depth_evidence
#   Step 03. base_depth_ratio
#   Step 04. plot_depth_ratio
#   Step 05. extract_sa_split_reads
#   Step 06. cluster_sa_split_reads
#   Step 07. extract_valid_tlen
#   Step 08. plot_tlen_distribution
#   Step 09. discordant_reads
#   Step 10. merge_evidence
#   Step 11. annotate_candidates
#
# Evidence merge logic:
#   split cluster is the candidate backbone;
#   depth and discordant evidence are used for validation.
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
# Directory and output path preparation
# ============================================================

sub prepare_sample_dirs {
    my ($outdir, $sample, $force) = @_;

    my %dir;

    $dir{sample}     = "$outdir/$sample";
    $dir{log}        = "$dir{sample}/00.log";
    $dir{target_bam} = "$dir{sample}/00.target_bam";
    $dir{depth}      = "$dir{sample}/01.depth";
    $dir{split}      = "$dir{sample}/02.split_reads";
    $dir{discordant} = "$dir{sample}/03.discordant_reads";
    $dir{candidate}  = "$dir{sample}/04.candidates";
    $dir{report}     = "$dir{sample}/05.report";
    $dir{tmp}        = "$dir{sample}/tmp";

    if (-d $dir{sample} && !$force) {
        die "[ERROR] Sample output directory already exists: $dir{sample}\n"
          . "        Please use --force to allow writing into an existing sample directory.\n";
    }

    for my $key (keys %dir) {
        make_path($dir{$key}) unless -d $dir{$key};
    }

    return %dir;
}

sub build_sample_outputs {
    my ($sample, $dir_ref) = @_;

    my $prefix = "$dir_ref->{depth}/$sample.depth_candidates";

    my %out = (
        target_bam             => "$dir_ref->{target_bam}/$sample.target.bam",
        target_bai             => "$dir_ref->{target_bam}/$sample.target.bam.bai",
        target_bed             => "$dir_ref->{target_bam}/$sample.target_regions.bed",
        target_region_tsv      => "$dir_ref->{target_bam}/$sample.target_regions.tsv",
        target_bam_summary     => "$dir_ref->{target_bam}/$sample.target_bam.summary.tsv",

        depth                  => "$prefix.tsv",
        all_window_ratio       => "$prefix.all_window_ratio.tsv",
        depth_plot_pdf         => "$dir_ref->{depth}/$sample.window_del.per_gene.pdf",
        gene_depth_dir         => "$prefix.gene_depth_files",

        cusum_dir              => "$prefix.cusum",
        cusum_all              => "$prefix.cusum.all.tsv",

        base_depth_ratio_dir   => "$prefix.base_depth_ratio",
        base_depth_ratio_all   => "$prefix.base_depth_ratio.all.tsv",
        base_depth_summary_all => "$prefix.base_depth.summary.all.tsv",

        split_raw              => "$dir_ref->{split}/$sample.split_reads.tsv",
        split_cluster          => "$dir_ref->{split}/$sample.split_reads.clusters.tsv",

        valid_tlen             => "$dir_ref->{discordant}/$sample.valid_tlen.tsv",
        tlen_plot_pdf          => "$dir_ref->{discordant}/$sample.tlen_distribution.pdf",
        discordant             => "$dir_ref->{discordant}/$sample.discordant_reads.tsv",

        merged                 => "$dir_ref->{candidate}/$sample.merged_candidates.tsv",
        annotated              => "$dir_ref->{report}/$sample.annotated_candidates.tsv",
    );

    return %out;
}

# ============================================================
# BAM list
# ============================================================

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
          . "        SampleID    /absolute/path/to/sample.bam\n"
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

        die "[ERROR] BAM path must be absolute at line $line_no:\n"
          . "        $bam\n"
            unless $bam =~ m{^/};

        die "[ERROR] BAM file does not end with .bam at line $line_no:\n"
          . "        $bam\n"
            unless $bam =~ /\.bam$/i;

        my $abs_bam = abs_path($bam);

        die "[ERROR] BAM file not found at line $line_no:\n"
          . "        $bam\n"
            unless defined $abs_bam && -s $abs_bam;

        check_bam_index($abs_bam);

        die "[ERROR] Duplicated sample name detected at line $line_no: $sample\n"
            if exists $seen_sample{$sample};

        die "[ERROR] Duplicated BAM path detected at line $line_no:\n"
          . "        $abs_bam\n"
            if exists $seen_bam{$abs_bam};

        $seen_sample{$sample} = 1;
        $seen_bam{$abs_bam}  = 1;

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
        unless -s $bai1 || -s $bai2;

    return 1;
}

# ============================================================
# Config and path helpers
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

        # Remove inline comments: KEY=value # comment
        $line =~ s/\s+#.*$//;

        next unless $line =~ /^\s*([^=\s]+)\s*=\s*(.*?)\s*$/;

        my $key = $1;
        my $val = $2;

        $key =~ s/^\s+|\s+$//g;
        $val =~ s/^\s+|\s+$//g;

        $val =~ s/^['"]//;
        $val =~ s/['"]$//;

        die "[ERROR] Empty config key in $file\n"
            if $key eq "";

        $conf{$key} = $val;
    }

    close $FH;

    return %conf;
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
    my ($v, $name) = @_;

    die "[ERROR] $name is undefined\n"
        unless defined $v;

    $v =~ s/^\s+|\s+$//g;

    return 1 if $v =~ /^(1|yes|true|on)$/i;
    return 0 if $v =~ /^(0|no|false|off)$/i;

    die "[ERROR] $name must be a boolean value: 1/0, yes/no, true/false, on/off. Observed: $v\n";
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

        die "[ERROR] Path not found: $path\n"
            unless defined $abs;

        return $abs;
    }

    my $full = "$project_root/$path";
    my $abs  = abs_path($full);

    die "[ERROR] Path not found: $full\n"
        unless defined $abs;

    return $abs;
}

# ============================================================
# Validation helpers
# ============================================================

sub validate_cusum_params {
    my ($cusum_ref) = @_;

    validate_baseline("CUSUM_BASELINE",              $cusum_ref->{baseline});
    validate_positive_integer("CUSUM_BIN_SIZE",       $cusum_ref->{bin_size});
    validate_nonnegative_number("CUSUM_K",            $cusum_ref->{k});
    validate_positive_number("CUSUM_H",               $cusum_ref->{h});
    validate_fraction("CUSUM_DEL_RATIO",              $cusum_ref->{del_ratio});
    validate_positive_integer("CUSUM_MIN_BINS",       $cusum_ref->{min_bins});
    validate_positive_integer("CUSUM_MIN_LEN",        $cusum_ref->{min_len});
    validate_fraction("CUSUM_EDGE_RATIO",             $cusum_ref->{edge_ratio});
    validate_fraction("CUSUM_RECOVER_RATIO",          $cusum_ref->{recover_ratio});
    validate_positive_integer("CUSUM_RECOVER_BINS",   $cusum_ref->{recover_bins});

    return 1;
}

sub validate_baseline {
    my ($name, $v) = @_;

    return 1 if defined $v && $v eq "auto";

    validate_positive_number($name, $v);

    return 1;
}

sub validate_fraction {
    my ($name, $v) = @_;

    die "[ERROR] $name must be numeric: $v\n"
        unless defined $v && $v =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;

    die "[ERROR] $name must be > 0 and <= 1: $v\n"
        unless $v > 0 && $v <= 1;

    return 1;
}

sub validate_positive_integer {
    my ($name, $v) = @_;

    die "[ERROR] $name must be a positive integer: $v\n"
        unless defined $v && $v =~ /^\d+$/ && $v >= 1;

    return 1;
}

sub validate_nonnegative_integer {
    my ($name, $v) = @_;

    die "[ERROR] $name must be a non-negative integer: $v\n"
        unless defined $v && $v =~ /^\d+$/ && $v >= 0;

    return 1;
}

sub validate_integer {
    my ($name, $v) = @_;

    die "[ERROR] $name must be an integer: $v\n"
        unless defined $v && $v =~ /^-?\d+$/;

    return 1;
}

sub validate_positive_number {
    my ($name, $v) = @_;

    die "[ERROR] $name must be numeric: $v\n"
        unless defined $v && $v =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;

    die "[ERROR] $name must be > 0: $v\n"
        unless $v > 0;

    return 1;
}

sub validate_nonnegative_number {
    my ($name, $v) = @_;

    die "[ERROR] $name must be numeric: $v\n"
        unless defined $v && $v =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;

    die "[ERROR] $name must be >= 0: $v\n"
        unless $v >= 0;

    return 1;
}

sub check_executable {
    my ($cmd, $name) = @_;

    die "[ERROR] Empty executable for $name\n"
        unless defined $cmd && $cmd ne "";

    if ($cmd =~ m{/}) {
        die "[ERROR] $name executable not found or not executable: $cmd\n"
            unless -x $cmd;
    }
    else {
        my $ret = system("command -v $cmd >/dev/null 2>&1");

        die "[ERROR] $name executable not found in PATH: $cmd\n"
            if $ret != 0;
    }

    return 1;
}

sub check_script {
    my ($path, $name) = @_;

    die "[ERROR] Pipeline script not found for $name: $path\n"
        unless defined $path && -s $path;

    return 1;
}

sub check_exon_txt_format {
    my ($file) = @_;

    my @required = qw/Gene Transcript Exon Chrom Start End Strand/;

    check_table_header(
        file     => $file,
        desc     => "REFSEQ_MANE_SELECT_EXON_TXT",
        required => \@required,
    );
}

sub check_gene_txt_format {
    my ($file) = @_;

    my @required = qw/Gene Transcript Chrom Start End Strand/;

    check_table_header(
        file     => $file,
        desc     => "REFSEQ_MANE_SELECT_GENE_TXT",
        required => \@required,
    );
}

sub check_table_header {
    my %args = @_;

    my $file     = $args{file};
    my $desc     = $args{desc};
    my $required = $args{required};

    open my $FH, "<", $file
        or die "[ERROR] Cannot open $desc file: $file\n";

    my $header = <$FH>;

    close $FH;

    die "[ERROR] Empty $desc file: $file\n"
        unless defined $header;

    chomp $header;
    $header =~ s/\r$//;

    my @cols = split /\t/, $header, -1;
    my %seen = map { $_ => 1 } @cols;

    for my $col (@$required) {
        die "[ERROR] Required column '$col' not found in $desc: $file\n"
            unless exists $seen{$col};
    }

    return 1;
}

sub check_core_gene_list_format {
    my ($file) = @_;

    open my $FH, "<", $file
        or die "[ERROR] Cannot open HCM_CORE_GENE_LIST: $file\n";

    my $line_no      = 0;
    my $data_line_no = 0;
    my $gene_count   = 0;

    while (my $line = <$FH>) {
        chomp $line;
        $line =~ s/\r$//;
        $line_no++;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        $line =~ s/\s+#.*$//;
        $line =~ s/^\s+|\s+$//g;

        next if $line eq "";

        $data_line_no++;

        my @f = split /\s+/, $line;
        my $gene = $f[0];
        $gene = "" unless defined $gene;

        if ($data_line_no == 1 && $gene =~ /^Gene$/i) {
            next;
        }

        die "[ERROR] Empty gene symbol in HCM_CORE_GENE_LIST at line $line_no\n"
            if $gene eq "";

        die "[ERROR] Invalid gene symbol in HCM_CORE_GENE_LIST at line $line_no: $gene\n"
            unless $gene =~ /^[A-Za-z0-9_.-]+$/;

        $gene_count++;
    }

    close $FH;

    die "[ERROR] No core genes found in HCM_CORE_GENE_LIST: $file\n"
        unless $gene_count > 0;

    return 1;
}

# ============================================================
# Shell and output helpers
# ============================================================

sub shell_quote {
    my ($str) = @_;

    die "[ERROR] Undefined shell argument\n"
        unless defined $str;

    return "''" if $str eq "";

    if ($str =~ /^[A-Za-z0-9_\.\-\/:=,+]+$/) {
        return $str;
    }

    $str =~ s/'/'\\''/g;

    return "'$str'";
}

sub print_summary {
    my %args = @_;

    print "\n";
    print "============================================================\n";
    print "HCMExonDel pipeline setup\n";
    print "============================================================\n";
    print "Config             : $args{config}\n";
    print "Project root       : $args{project_root}\n";
    print "BAM list           : $args{bam_list}\n";
    print "Exon TXT           : $args{exon_txt}\n";
    print "Gene TXT           : $args{gene_txt}\n";
    print "HCM core gene list : $args{core_gene_list}\n";
    print "Output directory   : $args{outdir}\n";
    print "PERL               : $args{perl}\n";
    print "PYTHON             : $args{python}\n";
    print "Sample number      : $args{sample_num}\n";

    print "\n";
    print "Workflow switches\n";
    print "------------------------------------------------------------\n";
    print "KEEP_GENE_DEPTH_FILE        : $args{keep_gene_depth_file}\n";
    print "OUTPUT_BASE_DEPTH_RATIO     : $args{output_base_depth_ratio}\n";
    print "OUTPUT_BASE_DEPTH_RATIO_ALL : $args{output_base_depth_ratio_all}\n";
    print "RUN_TLEN_QC                 : $args{run_tlen_qc}\n";
    print "KEEP_TMP                    : $args{keep_tmp}\n";

    print "\n";
    print "Target BAM and read extraction parameters\n";
    print "------------------------------------------------------------\n";
    print "TARGET_REGION_FLANK     : $args{target_region_flank}\n";
    print "TARGET_BAM_THREADS      : $args{target_bam_threads}\n";
    print "SPLIT_READ_THREADS      : $args{split_read_threads}\n";
    print "DISCORDANT_READ_THREADS : $args{discordant_read_threads}\n";

    if ($args{run_tlen_qc}) {
        my $t = $args{tlen_ref};

        print "\n";
        print "TLEN QC parameters\n";
        print "------------------------------------------------------------\n";
        print "TLEN_MAX_PAIRS   : $t->{max_pairs}\n";
        print "TLEN_THREADS     : $t->{threads}\n";
        print "TLEN_RANDOM_SEED : $t->{random_seed}\n";
    }

    print "\n";
    print "Evidence merging\n";
    print "------------------------------------------------------------\n";
    print "Merge backbone            : split-read cluster\n";
    print "Supporting evidence       : CUSUM depth + discordant read pairs\n";
    print "EVIDENCE_OVERLAP_FRACTION : $args{evidence_overlap_fraction}\n";

    print "\n";
    print "CUSUM parameters\n";
    print "------------------------------------------------------------\n";

    my $c = $args{cusum_ref};

    print "CUSUM_BASELINE      : $c->{baseline}\n";
    print "CUSUM_BIN_SIZE      : $c->{bin_size}\n";
    print "CUSUM_K             : $c->{k}\n";
    print "CUSUM_H             : $c->{h}\n";
    print "CUSUM_DEL_RATIO     : $c->{del_ratio}\n";
    print "CUSUM_MIN_BINS      : $c->{min_bins}\n";
    print "CUSUM_MIN_LEN       : $c->{min_len}\n";
    print "CUSUM_EDGE_RATIO    : $c->{edge_ratio}\n";
    print "CUSUM_RECOVER_RATIO : $c->{recover_ratio}\n";
    print "CUSUM_RECOVER_BINS  : $c->{recover_bins}\n";
    print "============================================================\n\n";
}

sub usage {
    return <<'USAGE';
Usage:
  perl HCMExonDel.pl \
    --config conf/hcm_exondel.example.conf \
    --bam-list bam.list \
    --outdir results

Required:
  --config      HCMExonDel config file.
  --bam-list    TAB-delimited BAM list with two columns:
                  SampleID    /absolute/path/to/sample.bam

Optional:
  --outdir      Output directory. Default: results
  --force       Allow writing into existing sample output directories.
  --help        Show this help message.

Generated output:
  OUTDIR/SAMPLE/SAMPLE.run.sh

Workflow:
  Step 00. extract_target_bam
  Step 01. gene_mean_depth
  Step 02. cusum_depth_evidence
  Step 03. base_depth_ratio
  Step 04. plot_depth_ratio
  Step 05. extract_sa_split_reads
  Step 06. cluster_sa_split_reads
  Step 07. extract_valid_tlen
  Step 08. plot_tlen_distribution
  Step 09. discordant_reads
  Step 10. merge_evidence
  Step 11. annotate_candidates

Notes:
  1. HCMExonDel only analyzes HCM core genes.
  2. The original WGS BAM is used only to generate SAMPLE.target.bam.
  3. All downstream analyses use SAMPLE.target.bam.
  4. merge_evidence.pl uses split clusters as candidate events and validates them
     with depth evidence and discordant read-pair evidence.

Example:
  perl HCMExonDel.pl \
    --config conf/hcm_exondel.example.conf \
    --bam-list test/bam.list \
    --outdir test/test_results \
    --force

Then run:
  bash test/test_results/SAMPLE/SAMPLE.run.sh
USAGE
}

