# HCMExonDel

HCMExonDel is a lightweight and configurable pipeline for detecting exon-level deletions in hypertrophic cardiomyopathy-related genes from whole-genome sequencing BAM files.

The pipeline integrates three complementary signals, including read-depth reduction, split-read evidence, and discordant read-pair evidence, to prioritize high-confidence exon-level deletion candidates for manual review and experimental validation.

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Workflow](#workflow)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input Files](#input-files)
- [Configuration](#configuration)
- [Quick Start](#quick-start)
- [Output Structure](#output-structure)
- [Output Files](#output-files)
- [Candidate Interpretation](#candidate-interpretation)
- [Manual Review and Validation](#manual-review-and-validation)
- [Running Individual Modules](#running-individual-modules)
- [Troubleshooting](#troubleshooting)
- [Limitations](#limitations)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

## Overview

Hypertrophic cardiomyopathy is a genetically heterogeneous cardiovascular disease. In some families, routine short-variant analysis or whole-exome sequencing may fail to identify disease-causing variants. Exon-level deletions and other structural variants can be missed by standard SNV/indel-oriented analysis pipelines.

HCMExonDel was developed to detect candidate exon-level deletions from WGS data, especially in HCM-related genes. It combines depth-ratio analysis with breakpoint-level evidence from split reads and discordant read pairs.

The pipeline is designed for research use and is suitable for WGS-based reanalysis of genetically unresolved HCM families.

## Key Features

- Detects exon-level deletion candidates from coordinate-sorted WGS BAM files.
- Supports both HCM core gene-focused analysis and all-gene analysis based on MANE RefSeq gene annotations.
- Uses three independent evidence types:
  - read-depth reduction,
  - SA-tag / split-read / soft-clipped read evidence,
  - discordant read-pair evidence.
- Prioritizes high-confidence candidates supported by all three evidence types.
- Generates an auxiliary all-evidence file for debugging and exploratory review.
- Annotates candidate deletions against exon-level MANE RefSeq annotations.
- Extracts candidate gene-level BAM files for IGV visualization.
- Uses a centralized configuration file.
- Generates one independent shell script for each sample.
- Does not require BED files; annotation intervals are provided as 1-based closed TXT files.

## Workflow

HCMExonDel uses a stepwise workflow.

```text
Input BAM
   |
   |-- Step 1: read-depth analysis
   |     run_gene_mean_depth.pl
   |     |
   |     |-- gene mean depth
   |     |-- window depth
   |     |-- depth ratio
   |     |-- deletion-supporting windows
   |
   |-- Step 1b: depth-ratio visualization
   |     plot_depth_ratio.py
   |
   |-- Step 2: split-read extraction
   |     extract_sa_split_reads.pl
   |
   |-- Step 3: split-read clustering
   |     cluster_sa_split_reads.pl
   |
   |-- Step 4: discordant read-pair analysis
   |     run_discordant_reads.pl
   |
   |-- Step 5: evidence merging
   |     merge_evidence.pl
   |
   |-- Step 6: candidate gene BAM extraction
   |     extract_candidate_gene_bam.pl
   |
   |-- Step 7: candidate annotation
         annotate_candidates.pl
```

The main script `HCMExonDel.pl` does not directly execute all analysis steps. Instead, it generates a sample-specific shell script:

```text
OUTDIR/SAMPLE/SAMPLE.run.sh
```

The generated shell script can then be executed manually or submitted to an HPC scheduler according to the user's environment.

## Evidence Types

### 1. Read-depth evidence

Read-depth evidence is based on reduced sequencing depth across consecutive genomic windows within a gene region.

The pipeline calculates:

- base-level depth,
- gene-level mean depth,
- window-level mean depth,
- window-level depth ratio,
- consecutive deletion-supporting windows.

A depth candidate is generated when multiple consecutive windows show a depth ratio below the configured deletion cutoff.

### 2. Split-read evidence

Split-read evidence is extracted from reads with supplementary alignment information or soft-clipped sequence signals. These reads may support deletion breakpoints and can help refine the candidate interval.

The pipeline first extracts raw split-read signals and then clusters them into candidate deletion-supporting events.

### 3. Discordant read-pair evidence

Discordant read-pair evidence is based on paired-end reads with abnormal insert size or deletion-like orientation. These read pairs can support structural variation spanning the candidate interval.

### 4. Evidence integration

`merge_evidence.pl` integrates the three evidence types by gene and chromosome.

The main merged candidate output keeps only candidates supported by all three evidence types:

```text
Depth + Split + Discordant
```

An additional `.all.tsv` file is also generated, which contains all evidence-supported regions with at least one evidence type. This file is useful for debugging, threshold tuning, and exploratory review.

## Repository Structure

```text
HCMExonDel/
├── HCMExonDel.pl
├── bin/
│   ├── annotate_candidates.pl
│   ├── cluster_sa_split_reads.pl
│   ├── extract_candidate_gene_bam.pl
│   ├── extract_sa_split_reads.pl
│   ├── merge_evidence.pl
│   ├── plot_depth_ratio.py
│   ├── run_discordant_reads.pl
│   ├── run_gene_mean_depth.pl
│   ├── run_split_reads.pl
│   └── simulate_gene_wgs_del.pl
├── conf/
│   └── hcm_exondel.example.conf
├── db/
│   ├── canonical_transcripts/
│   ├── exon_annotation_bed/
│   ├── hcm_core_genes.txt
│   └── README
├── example/
├── test/
└── .gitignore
```

### Main script

| File | Description |
|---|---|
| `HCMExonDel.pl` | Main driver script. It checks inputs and generates one `SAMPLE.run.sh` file for each BAM sample. |

### Core scripts

| File | Description |
|---|---|
| `bin/run_gene_mean_depth.pl` | Performs gene-level and window-level depth-ratio analysis. |
| `bin/plot_depth_ratio.py` | Generates per-gene depth-ratio PDF plots. |
| `bin/extract_sa_split_reads.pl` | Extracts SA-tag, split-read, and soft-clipped read signals. |
| `bin/cluster_sa_split_reads.pl` | Clusters split-read signals into candidate breakpoint events. |
| `bin/run_discordant_reads.pl` | Detects and clusters discordant read-pair signals. |
| `bin/merge_evidence.pl` | Merges depth, split-read, and discordant-read evidence. |
| `bin/extract_candidate_gene_bam.pl` | Extracts gene-level BAM files for candidate genes. |
| `bin/annotate_candidates.pl` | Annotates merged candidates against exon annotation. |
| `bin/simulate_gene_wgs_del.pl` | Utility script for simulating gene-level WGS deletion signals. |

## Requirements

### Operating system

HCMExonDel is designed for Linux-based systems.

It can be used on:

- local Linux workstations,
- Linux servers,
- HPC clusters.

### Required software

- Perl
- Python 3
- samtools

### Recommended software

- IGV
- bcftools
- GNU coreutils
- standard Linux shell environment

### Perl modules

Most scripts use standard Perl modules, including:

```text
Getopt::Long
File::Basename
File::Path
Cwd
FindBin
```

### Python packages

The plotting script may require:

```text
pandas
numpy
matplotlib
```

Install them if needed:

```bash
pip install pandas numpy matplotlib
```

## Installation

Clone the repository:

```bash
git clone https://github.com/Xiaohuaniu0032/HCMExonDel.git
cd HCMExonDel
```

Check the main script:

```bash
perl HCMExonDel.pl --help
```

Check samtools:

```bash
samtools --version
```

Optional: make scripts executable.

```bash
chmod +x HCMExonDel.pl
chmod +x bin/*.pl
chmod +x bin/*.py
```

## Input Files

### 1. BAM file

Each sample should have a coordinate-sorted WGS BAM file.

Requirements:

- BAM file must be coordinate-sorted.
- BAM file must have an index file.
- BAM path must be an absolute path.
- BAM filename must end with `.bam`.
- Chromosome naming must be consistent with the annotation files.

Valid BAM index formats:

```text
sample.bam.bai
sample.bai
```

### 2. BAM list file

The BAM list file must contain two TAB-delimited columns.

```text
SampleName    /absolute/path/to/sample.bam
```

Example:

```text
25B09089386    /data/project/HCM/25B09089386.final.merge.bam
25B09089387    /data/project/HCM/25B09089387.final.merge.bam
```

Rules:

- No header is required.
- Empty lines are ignored.
- Lines beginning with `#` are ignored.
- Sample names can contain letters, numbers, dots, underscores, and hyphens.
- BAM paths must be absolute paths.
- Duplicate sample names are not allowed.
- Duplicate BAM paths are not allowed.

### 3. MANE RefSeq exon annotation file

The exon annotation file is configured by:

```text
REFSEQ_MANE_SELECT_EXON_TXT
```

Required columns:

```text
Gene
Transcript
Exon
Chrom
Start
End
Strand
```

Coordinates are 1-based closed intervals.

### 4. MANE RefSeq gene annotation file

The gene annotation file is configured by:

```text
REFSEQ_MANE_SELECT_GENE_TXT
```

Required columns:

```text
Gene
Transcript
Chrom
Start
End
Strand
ExonCount
```

Coordinates are 1-based closed intervals.

### 5. HCM core gene list

The HCM core gene list is configured by:

```text
HCM_CORE_GENE_LIST
```

The file should contain one gene symbol per line.

An optional second column can be used to store gene classification.

Example:

```text
MYBPC3    Definitive
MYH7      Definitive
FHOD3     Moderate
```

When the following option is enabled, only genes in this list are analyzed:

```text
ANALYZE_CORE_GENES_ONLY=1
```

When it is disabled, all genes in the MANE RefSeq gene annotation file are analyzed:

```text
ANALYZE_CORE_GENES_ONLY=0
```

## Configuration

The example configuration file is:

```text
conf/hcm_exondel.example.conf
```

Create your own configuration file before running:

```bash
cp conf/hcm_exondel.example.conf conf/hcm_exondel.conf
```

Then edit the software paths and parameter values.

### Important configuration parameters

#### Software

```text
PERL=/usr/bin/perl
SAMTOOLS=/path/to/samtools
```

#### Annotation files

```text
REFSEQ_MANE_SELECT_EXON_TXT=db/exon_annotation_bed/RefSeq_MANE_Select.exon.txt
REFSEQ_MANE_SELECT_GENE_TXT=db/exon_annotation_bed/RefSeq_MANE_Select.gene.txt
REF_FASTA_INDEX=
```

`REF_FASTA_INDEX` is optional but recommended for discordant read-pair analysis. If provided, genomic intervals can be clipped by chromosome length.

#### Target gene selection

```text
HCM_CORE_GENE_LIST=db/hcm_core_genes.txt
ANALYZE_CORE_GENES_ONLY=1
TARGET_REGION_FLANK=5000
```

`TARGET_REGION_FLANK` adds flanking sequence to target gene regions when scanning split reads and discordant read pairs. This is useful for capturing breakpoint-supporting reads near gene boundaries.

#### Common read filters

```text
MIN_MAPQ=20
EXCLUDE_DUPLICATES=1
```

`MIN_MAPQ` filters low mapping-quality reads.

`EXCLUDE_DUPLICATES=1` removes duplicate reads from read-depth, split-read, and discordant-read analysis.

#### Depth-ratio analysis

```text
WINDOW_SIZE=1000
WINDOW_STEP=500
MIN_WINDOW_SIZE=100
MIN_GENE_MEAN_DEPTH=20
DEL_DEPTH_RATIO_CUTOFF=0.65
MIN_CONSECUTIVE_DEL_WINDOWS=3
KEEP_GENE_DEPTH_FILE=1
```

`DEL_DEPTH_RATIO_CUTOFF=0.65` is a practical default for heterozygous exon-level deletions.

`KEEP_GENE_DEPTH_FILE=1` keeps and reuses per-gene base-level depth files.

#### Split-read clustering

```text
SA_SPLIT_CLUSTER_WINDOW=20
SA_SPLIT_MIN_SUPPORT_READS=5
```

`SA_SPLIT_CLUSTER_WINDOW` defines the maximum allowed breakpoint difference for clustering split-read signals.

`SA_SPLIT_MIN_SUPPORT_READS` defines the minimum number of unique read names required for a split-read cluster to pass.

#### Discordant read-pair analysis

```text
MIN_DISCORDANT_INSERT_SIZE=1000
MIN_DISCORDANT_READS=3
DISCORDANT_CLUSTER_DISTANCE=500
FILTER_DELETION_ORIENTATION=1
```

`FILTER_DELETION_ORIENTATION=1` keeps deletion-like inward-facing paired-end orientation.

#### Exon annotation

```text
MIN_EXON_OVERLAP_FRACTION=0.20
```

This parameter controls the minimum fraction of an exon that must be overlapped by a candidate interval to be considered affected.

#### Runtime settings

```text
KEEP_TMP=0
VERBOSE=1
```

`THREADS` is not used by the current pipeline.

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/Xiaohuaniu0032/HCMExonDel.git
cd HCMExonDel
```

### 2. Prepare the configuration file

```bash
cp conf/hcm_exondel.example.conf conf/hcm_exondel.conf
```

Edit the configuration file:

```bash
vim conf/hcm_exondel.conf
```

At minimum, check the following parameters:

```text
PERL
SAMTOOLS
REFSEQ_MANE_SELECT_EXON_TXT
REFSEQ_MANE_SELECT_GENE_TXT
HCM_CORE_GENE_LIST
```

### 3. Prepare the BAM list

Create a file named `input_bam.list`.

```text
25B09089386    /absolute/path/to/25B09089386.final.merge.bam
```

The file must be TAB-delimited.

### 4. Generate sample-specific shell scripts

```bash
perl HCMExonDel.pl \
  --config conf/hcm_exondel.conf \
  --bam-list input_bam.list \
  --outdir results
```

This command generates:

```text
results/25B09089386/25B09089386.run.sh
```

### 5. Run the sample shell script

```bash
bash results/25B09089386/25B09089386.run.sh
```

For multiple samples, run each generated shell script manually or submit them to your HPC scheduler.

### 6. Re-run into an existing sample directory

By default, the main script stops if the sample output directory already exists.

Use `--force` to allow writing into an existing sample directory:

```bash
perl HCMExonDel.pl \
  --config conf/hcm_exondel.conf \
  --bam-list input_bam.list \
  --outdir results \
  --force
```

## Output Structure

For each sample, HCMExonDel generates the following directory structure:

```text
OUTDIR/
└── SAMPLE/
    ├── SAMPLE.run.sh
    ├── 00.log/
    ├── 01.depth/
    ├── 02.split_reads/
    ├── 03.discordant_reads/
    ├── 04.candidates/
    ├── 05.gene_bam/
    ├── 06.report/
    └── tmp/
```

No `sample_shell.list` is generated.

No `qsub_command.list` is generated.

## Output Files

### 1. Log files

Directory:

```text
00.log/
```

Typical files:

```text
01.gene_mean_depth.log
01b.plot_depth_ratio.log
02.extract_sa_split_reads.log
03.cluster_sa_split_reads.log
04.discordant_reads.log
05.merge_evidence.log
06.extract_candidate_gene_bam.log
07.annotate_candidates.log
```

### 2. Depth analysis output

Directory:

```text
01.depth/
```

Typical files:

```text
SAMPLE.depth_candidates.tsv
SAMPLE.depth_candidates.all_window_ratio.tsv
SAMPLE.depth_candidates.del_windows.tsv
SAMPLE.depth_candidates.gene_mean_depth.tsv
SAMPLE.depth_candidates.window_depth.tsv
SAMPLE.window_del.per_gene.pdf
SAMPLE.depth_candidates.gene_depth_files/
```

Description:

| File | Description |
|---|---|
| `SAMPLE.depth_candidates.tsv` | Merged depth-supported deletion candidates. |
| `SAMPLE.depth_candidates.all_window_ratio.tsv` | All analyzed windows with depth ratios. |
| `SAMPLE.depth_candidates.del_windows.tsv` | Windows passing the deletion depth-ratio cutoff. |
| `SAMPLE.depth_candidates.gene_mean_depth.tsv` | Gene-level mean depth. |
| `SAMPLE.depth_candidates.window_depth.tsv` | Window-level depth summary. |
| `SAMPLE.window_del.per_gene.pdf` | Per-gene depth-ratio visualization. |
| `SAMPLE.depth_candidates.gene_depth_files/` | Per-gene base-level depth files. |

### 3. Split-read output

Directory:

```text
02.split_reads/
```

Typical files:

```text
SAMPLE.split_reads.tsv
SAMPLE.split_reads.clusters.tsv
SAMPLE.split_reads.failed_clusters.tsv
SAMPLE.split_reads.supporting_reads.tsv
```

Description:

| File | Description |
|---|---|
| `SAMPLE.split_reads.tsv` | Raw SA-tag, split-read, or soft-clipped read signals. |
| `SAMPLE.split_reads.clusters.tsv` | Passed split-read clusters. Used by `merge_evidence.pl`. |
| `SAMPLE.split_reads.failed_clusters.tsv` | Split-read clusters that did not pass filtering. |
| `SAMPLE.split_reads.supporting_reads.tsv` | Supporting read-level information for passed clusters. |

### 4. Discordant read-pair output

Directory:

```text
03.discordant_reads/
```

Typical files:

```text
SAMPLE.discordant_reads.tsv
SAMPLE.discordant_reads.raw_discordant_pairs.tsv
SAMPLE.discordant_reads.supporting_pairs.tsv
SAMPLE.discordant_reads.discarded_clusters.tsv
```

Description:

| File | Description |
|---|---|
| `SAMPLE.discordant_reads.tsv` | Passed discordant read-pair clusters. Used by `merge_evidence.pl`. |
| `SAMPLE.discordant_reads.raw_discordant_pairs.tsv` | Raw discordant read-pair records. |
| `SAMPLE.discordant_reads.supporting_pairs.tsv` | Supporting read-pair records for passed clusters. |
| `SAMPLE.discordant_reads.discarded_clusters.tsv` | Discordant clusters that did not pass filtering. |

### 5. Candidate merging output

Directory:

```text
04.candidates/
```

Typical files:

```text
SAMPLE.merged_candidates.tsv
SAMPLE.merged_candidates.tsv.all.tsv
```

Description:

| File | Description |
|---|---|
| `SAMPLE.merged_candidates.tsv` | Main candidate file. Keeps only candidates supported by depth, split-read, and discordant-read evidence. |
| `SAMPLE.merged_candidates.tsv.all.tsv` | Auxiliary file containing all evidence-supported regions with at least one evidence type. |

### 6. Candidate gene BAM output

Directory:

```text
05.gene_bam/
```

Typical files:

```text
GENE.bam
GENE.bam.bai
```

Candidate gene BAM files are extracted for genes appearing in the merged candidate file.

The current main script always generates this step. No `--flank` parameter is passed to `extract_candidate_gene_bam.pl`; therefore, gene BAM extraction uses the script's internal default flank value of `0` and is strictly based on gene coordinates.

### 7. Final annotation report

Directory:

```text
06.report/
```

Typical file:

```text
SAMPLE.annotated_candidates.tsv
```

This file contains exon-level annotation for merged candidate deletions.

## Candidate Interpretation

### Main candidate file

The main merged candidate file contains high-confidence candidates supported by all three evidence types.

Typical important columns include:

| Column | Description |
|---|---|
| `SampleID` | Sample ID. |
| `Gene` | Candidate gene. |
| `Classification` | Gene classification, if provided by the HCM gene list. |
| `Transcript` | Transcript ID. |
| `Chrom` | Chromosome. |
| `Merged_Start` | Outer start boundary of all overlapping evidence intervals. |
| `Merged_End` | Outer end boundary of all overlapping evidence intervals. |
| `Merged_Size` | Size of the merged outer interval. |
| `Core_Start` | Conservative evidence-overlap start position. |
| `Core_End` | Conservative evidence-overlap end position. |
| `Core_Size` | Size of the conservative evidence-overlap interval. |
| `Best_Start` | Recommended candidate start position. |
| `Best_End` | Recommended candidate end position. |
| `Best_Size` | Recommended candidate size. |
| `Best_Evidence` | Evidence type used to define the recommended interval. |
| `Evidence_Count` | Number of supporting evidence types. |
| `Evidence_Level` | Evidence confidence level. |
| `Evidence_Types` | Supporting evidence types. |
| `Depth_Support` | Whether depth evidence supports the candidate. |
| `Split_Read_Support` | Whether split-read evidence supports the candidate. |
| `Discordant_Read_Support` | Whether discordant read-pair evidence supports the candidate. |
| `Candidate_Status` | Candidate status based on evidence count. |

### Coordinate definitions

HCMExonDel reports three coordinate systems.

#### 1. `Merged_Start` and `Merged_End`

These represent the outer boundary of all overlapping evidence intervals.

This interval may be wider than the true deletion because it includes the full span of all supporting evidence.

#### 2. `Core_Start` and `Core_End`

These represent the conservative genomic segment supported by the listed evidence types.

For three-evidence candidates, the core interval is the intersection of depth, split-read, and discordant-read evidence.

#### 3. `Best_Start` and `Best_End`

These are the recommended candidate boundaries for downstream annotation and reporting.

Boundary priority:

```text
Split-read evidence > Discordant read-pair evidence > Depth evidence
```

Rationale:

- Split-read clusters are usually closest to the real breakpoint.
- Discordant read-pair clusters provide approximate breakpoint support.
- Depth intervals are useful but usually have coarser boundaries because of window-based analysis.

## Manual Review and Validation

Candidate deletions should be manually reviewed before interpretation.

Recommended review steps:

1. Load the original BAM file in IGV.
2. Load the candidate gene-level BAM file from `05.gene_bam/`.
3. Inspect read depth across the candidate region.
4. Check whether the target exon region shows a consistent depth reduction.
5. Inspect soft-clipped or split-read signals near predicted breakpoints.
6. Inspect discordant read pairs spanning the candidate interval.
7. Compare with family members or control samples if available.
8. Prioritize candidates affecting coding exons and segregating with disease.

Recommended validation methods include:

- ddPCR,
- MLPA,
- breakpoint PCR,
- Sanger sequencing across the breakpoint,
- RT-PCR or RNA-seq for transcript-level confirmation,
- family segregation analysis.

## Running Individual Modules

The recommended way is to run the main script first and then execute the generated sample shell script.

However, each module can also be run individually for debugging.

### Depth analysis

```bash
perl bin/run_gene_mean_depth.pl \
  --config conf/hcm_exondel.conf \
  --bam /absolute/path/to/sample.bam \
  --sample SAMPLE \
  --out results/SAMPLE/01.depth/SAMPLE.depth_candidates.tsv
```

### Depth-ratio plotting

```bash
python3 bin/plot_depth_ratio.py \
  --conf conf/hcm_exondel.conf \
  --input results/SAMPLE/01.depth/SAMPLE.depth_candidates.all_window_ratio.tsv \
  --output results/SAMPLE/01.depth/SAMPLE.window_del.per_gene.pdf
```

### Split-read extraction

```bash
perl bin/extract_sa_split_reads.pl \
  --conf conf/hcm_exondel.conf \
  --bam /absolute/path/to/sample.bam \
  --sample SAMPLE \
  --out results/SAMPLE/02.split_reads/SAMPLE.split_reads.tsv
```

### Split-read clustering

```bash
perl bin/cluster_sa_split_reads.pl \
  --conf conf/hcm_exondel.conf \
  --input results/SAMPLE/02.split_reads/SAMPLE.split_reads.tsv \
  --outfile results/SAMPLE/02.split_reads/SAMPLE.split_reads.clusters.tsv
```

### Discordant read-pair analysis

```bash
perl bin/run_discordant_reads.pl \
  --config conf/hcm_exondel.conf \
  --bam /absolute/path/to/sample.bam \
  --sample SAMPLE \
  --out results/SAMPLE/03.discordant_reads/SAMPLE.discordant_reads.tsv
```

### Evidence merging

```bash
perl bin/merge_evidence.pl \
  --config conf/hcm_exondel.conf \
  --sample SAMPLE \
  --depth results/SAMPLE/01.depth/SAMPLE.depth_candidates.tsv \
  --split results/SAMPLE/02.split_reads/SAMPLE.split_reads.clusters.tsv \
  --discordant results/SAMPLE/03.discordant_reads/SAMPLE.discordant_reads.tsv \
  --out results/SAMPLE/04.candidates/SAMPLE.merged_candidates.tsv
```

Important: `merge_evidence.pl` requires `SAMPLE.split_reads.clusters.tsv`, not the raw `SAMPLE.split_reads.tsv` file.

### Candidate gene BAM extraction

```bash
perl bin/extract_candidate_gene_bam.pl \
  --config conf/hcm_exondel.conf \
  --candidate results/SAMPLE/04.candidates/SAMPLE.merged_candidates.tsv \
  --bam /absolute/path/to/sample.bam \
  --outdir results/SAMPLE/05.gene_bam
```

### Candidate annotation

```bash
perl bin/annotate_candidates.pl \
  --config conf/hcm_exondel.conf \
  --input results/SAMPLE/04.candidates/SAMPLE.merged_candidates.tsv \
  --out results/SAMPLE/06.report/SAMPLE.annotated_candidates.tsv
```

## Troubleshooting

### 1. BAM index not found

Error example:

```text
[ERROR] BAM index not found
```

Generate the BAM index:

```bash
samtools index sample.bam
```

Accepted index files:

```text
sample.bam.bai
sample.bai
```

### 2. Invalid BAM list format

The BAM list must contain exactly two TAB-delimited columns.

Correct:

```text
SampleA    /absolute/path/to/SampleA.bam
```

Incorrect:

```text
SampleA /absolute/path/to/SampleA.bam
```

Make sure the separator is a TAB, not spaces.

### 3. Output directory already exists

By default, HCMExonDel stops when the sample output directory already exists.

Use:

```bash
--force
```

Example:

```bash
perl HCMExonDel.pl \
  --config conf/hcm_exondel.conf \
  --bam-list input_bam.list \
  --outdir results \
  --force
```

### 4. Chromosome naming mismatch

If BAM uses `chr18` but annotation files use `18`, or vice versa, region extraction may fail.

Make sure the BAM, reference genome, FASTA index, gene annotation file, and exon annotation file use the same chromosome naming convention.

### 5. No depth candidate detected

Possible reasons:

- No exon-level deletion exists in the analyzed genes.
- Sequencing depth is too low or uneven.
- `DEL_DEPTH_RATIO_CUTOFF` is too strict.
- `MIN_CONSECUTIVE_DEL_WINDOWS` is too high.
- `WINDOW_SIZE` and `WINDOW_STEP` are not suitable for the exon size.
- The gene is not included when `ANALYZE_CORE_GENES_ONLY=1`.

### 6. No split-read cluster detected

Possible reasons:

- The true breakpoint is in a repetitive or poorly mappable region.
- SA tags are absent or filtered during BAM processing.
- `SA_SPLIT_MIN_SUPPORT_READS` is too strict.
- `MIN_MAPQ` is too strict.
- The deletion is supported by depth but not breakpoint-spanning reads.

### 7. Weak discordant read-pair evidence

Possible reasons:

- The deletion size is smaller than the insert-size threshold.
- `MIN_DISCORDANT_INSERT_SIZE` is too high.
- `MIN_DISCORDANT_READS` is too strict.
- Library insert-size distribution is different from the default assumption.
- Local mapping ambiguity reduces paired-end support.

### 8. Merged candidate file is empty

The main merged candidate file keeps only candidates supported by all three evidence types.

Check the auxiliary file:

```text
SAMPLE.merged_candidates.tsv.all.tsv
```

This file contains all regions with at least one evidence type and can help determine whether thresholds are too strict.

## Limitations

- HCMExonDel is designed for WGS data and is not optimized for WES data.
- Depth-ratio analysis depends on sequencing depth, GC bias, mappability, and local alignment quality.
- Split-read and discordant read-pair signals may be weak or absent in repetitive regions.
- Small deletions may be difficult to detect using window-based depth analysis.
- The main candidate file is intentionally stringent and keeps only three-evidence candidates.
- The pipeline does not perform clinical variant classification.
- All candidates require manual review and independent experimental validation.

## Recommended Practices

- Use high-quality, coordinate-sorted WGS BAM files.
- Use annotation files matching the same reference genome as the BAM.
- Keep chromosome naming consistent across all files.
- Start with the default configuration and then tune thresholds based on sequencing depth and library characteristics.
- Review both the main candidate file and the `.all.tsv` auxiliary file.
- Validate high-confidence candidates using orthogonal experimental methods.
- For family studies, evaluate segregation between candidate deletions and disease status.

## Citation

If you use HCMExonDel in a publication, report, or internal research project, please cite this repository and describe the main evidence-integration strategy.

Suggested method description:

```text
Exon-level deletion candidates in hypertrophic cardiomyopathy-related genes were detected from WGS BAM files using HCMExonDel, a custom pipeline integrating read-depth reduction, SA-tag/split-read evidence, and discordant read-pair evidence. Candidate intervals supported by multiple evidence types were annotated against MANE RefSeq exon annotations and manually reviewed using IGV.
```

## License

No license file is currently provided.

Please add a license according to the intended use of this repository, for example:

```text
MIT License
Apache License 2.0
GPL-3.0
Research-use-only license
```

## Contact

Maintainer:

```text
Xiaohuaniu0032
GitHub: https://github.com/Xiaohuaniu0032/HCMExonDel
```

## Disclaimer

HCMExonDel is provided for research use only. It is not intended for direct clinical diagnosis. Candidate deletion events should be reviewed by qualified personnel and validated using independent experimental methods before clinical interpretation.
