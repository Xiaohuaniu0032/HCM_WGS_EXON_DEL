# HCMExonDel

HCMExonDel is a lightweight, configurable pipeline for detecting candidate exon-level deletions in hypertrophic cardiomyopathy (HCM)-related core genes from whole-genome sequencing (WGS) BAM files.

The current workflow is **split-read-centered**: split-read clusters define candidate deletion intervals, while CUSUM depth evidence and discordant read-pair evidence are used to validate and prioritize those intervals. Final candidates are annotated against MANE RefSeq exon coordinates for manual review and experimental validation.

> **Research use only.** HCMExonDel does not perform clinical variant classification, and all reported candidates require independent review and validation.

---

## Table of Contents

- [Overview](#overview)
- [Current Design](#current-design)
- [Workflow](#workflow)
- [Evidence Integration](#evidence-integration)
- [Requirements](#requirements)
- [Repository Structure](#repository-structure)
- [Input Files](#input-files)
- [Configuration](#configuration)
- [Quick Start](#quick-start)
- [Output Structure](#output-structure)
- [Key Output Files](#key-output-files)
- [Candidate Interpretation](#candidate-interpretation)
- [Test Example: 25B09089386](#test-example-25b09089386)
- [Manual Review and Validation](#manual-review-and-validation)
- [Running Individual Modules](#running-individual-modules)
- [Troubleshooting](#troubleshooting)
- [Limitations](#limitations)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## Overview

Exon-level deletions and other structural variants may be missed by analysis pipelines focused mainly on SNVs and small indels. HCMExonDel was developed for WGS-based reanalysis of genetically unresolved HCM cases and families, with emphasis on deletions affecting established HCM-related genes.

HCMExonDel integrates three complementary signals:

1. **Read-depth reduction** detected by base-level CUSUM analysis.
2. **Split-read evidence** derived from supplementary alignments and breakpoint-supporting reads.
3. **Discordant read-pair evidence** derived from abnormally large, deletion-like paired-end fragments.

The pipeline does not call deletions across the whole genome. It always analyzes genes listed in `HCM_CORE_GENE_LIST`.

---

## Current Design

The current implementation differs from earlier versions of the pipeline in several important ways.

### Core-gene-only analysis

HCMExonDel always analyzes the HCM core genes listed in:

```text
HCM_CORE_GENE_LIST
```

The former parameter `ANALYZE_CORE_GENES_ONLY` has been removed. If it is still present in a configuration file, the main script stops with an error and asks the user to remove it.

### Target BAM first

For each sample, the original WGS BAM is used only in Step 00 to generate:

```text
SAMPLE.target.bam
```

The target BAM contains the core-gene regions plus `TARGET_REGION_FLANK`. All downstream depth, split-read, TLEN, and discordant-read analyses use this smaller target BAM.

### Split-read-centered candidate generation

A passed split-read cluster is the candidate event backbone. CUSUM depth intervals and discordant-read clusters are evaluated as supporting evidence for that split interval.

Consequently, a depth-only or discordant-only signal is not included in the final merged candidate file when no qualifying split-read cluster is present.

### Candidate gene BAM extraction removed

Candidate gene BAM extraction is no longer part of the standard workflow. The former parameter `CANDIDATE_BAM_FLANK` and the former `05.gene_bam/` output directory have been removed.

For IGV review, use the original WGS BAM and/or the generated `SAMPLE.target.bam`.

---

## Workflow

The main script, `HCMExonDel.pl`, validates the configuration and input BAM list, then generates one executable shell script for each sample. It does not directly execute the complete analysis.

```text
Original coordinate-sorted WGS BAM
                |
                v
Step 00  extract_target_bam
         Core genes + flanking regions -> SAMPLE.target.bam
                |
                v
Step 01  gene_mean_depth
         Gene depth, window depth, depth ratios, low-depth windows
                |
                v
Step 02  cusum_depth_evidence
         Base-level CUSUM deletion-supporting intervals
                |
                +------------------------------+
                |                              |
                v                              v
Step 03  base_depth_ratio              Step 04  plot_depth_ratio
         Optional diagnostic output             Per-gene PDF
                |
                v
Step 05  extract_sa_split_reads
         Raw breakpoint-supporting records
                |
                v
Step 06  cluster_sa_split_reads
         Passed split-read clusters
                |
                +------------------------------+
                |                              |
                v                              v
Step 07  extract_valid_tlen             Step 09  discordant_reads
Step 08  plot_tlen_distribution                 Passed discordant clusters
         Optional library QC
                \                              /
                 \                            /
                  v                          v
Step 10             merge_evidence
        Split clusters validated by CUSUM depth and discordant evidence
                                |
                                v
Step 11             annotate_candidates
        MANE RefSeq exon annotation and final report
```

For each sample, the main script generates:

```text
OUTDIR/SAMPLE/SAMPLE.run.sh
```

Run this shell script manually or submit it to an HPC scheduler.

---

## Evidence Integration

### 1. Window-level depth analysis

`run_gene_mean_depth.pl` calculates gene-level and window-level depth statistics using the target BAM. A window is considered deletion-supporting when its depth ratio is below `DEL_DEPTH_RATIO_CUTOFF`.

The window-level results are useful for:

- detecting consecutive low-depth windows;
- visualizing depth changes across each core gene;
- checking whether a candidate interval shows an expected heterozygous depth reduction.

The window candidate file is **not** the depth input used by the current evidence-merging step.

### 2. CUSUM depth evidence

`run_sample_cusum_depth.pl` applies `cusum_depth_del.pl` to retained per-gene base-depth files. The combined output is:

```text
SAMPLE.depth_candidates.cusum.all.tsv
```

This file is the depth-evidence input passed to `merge_evidence.pl`.

For each split interval, overlapping CUSUM depth intervals are combined. Depth support is assigned when their covered portion reaches `EVIDENCE_OVERLAP_FRACTION` of the split interval.

### 3. Base-depth ratio output

When `OUTPUT_BASE_DEPTH_RATIO=1`, the pipeline additionally calculates base-level depth-ratio files for inspection and debugging.

These outputs are diagnostic and are not used by the current `merge_evidence.pl` step.

### 4. Split-read evidence

`extract_sa_split_reads.pl` extracts breakpoint-supporting records from the target BAM, and `cluster_sa_split_reads.pl` groups nearby signals into candidate events.

A cluster must contain at least `SA_SPLIT_MIN_SUPPORT_READS` unique supporting read names to pass.

The passed cluster file is:

```text
SAMPLE.split_reads.clusters.tsv
```

This file defines the candidate intervals evaluated in Step 10.

### 5. Discordant read-pair evidence

`run_discordant_reads.pl` identifies paired-end fragments with an insert size at least `MIN_DISCORDANT_INSERT_SIZE`, optionally requiring deletion-like orientation. Nearby pairs are clustered, and a cluster must contain at least `MIN_DISCORDANT_READS` supporting read pairs to pass.

A split candidate receives discordant support when a qualifying overlapping discordant cluster reaches the configured overlap requirement.

### 6. Evidence levels

The current confidence model is:

| Evidence level | Required evidence | Interpretation |
|---|---|---|
| `High` | Split + Depth + Discordant | Breakpoint, depth, and paired-end evidence are all present. |
| `Moderate` | Split + Depth, without qualifying Discordant support | Breakpoint and depth evidence are present. |
| `Low` | Split without qualifying Depth support | A split cluster is present, but depth support does not pass the overlap threshold. Discordant evidence alone does not raise the candidate to `Moderate`. |

The final merged file contains all passed split-centered candidates and records their evidence level. The current workflow does not generate a separate `.all.tsv` candidate file.

---

## Requirements

### Operating system

HCMExonDel is designed for Linux environments, including Linux workstations, servers, and HPC clusters.

### Required software

- Perl 5
- Python 3
- samtools
- Bash

### Perl modules

The pipeline mainly uses standard Perl modules, including:

```text
Getopt::Long
File::Basename
File::Path
Cwd
FindBin
```

### Python packages

The plotting scripts require:

```text
pandas
numpy
matplotlib
```

Install them when necessary:

```bash
python3 -m pip install pandas numpy matplotlib
```

### Recommended review software

- IGV

---

## Repository Structure

```text
HCMExonDel/
├── HCMExonDel.pl
├── bin/
│   ├── extract_target_bam.pl
│   ├── run_gene_mean_depth.pl
│   ├── run_sample_cusum_depth.pl
│   ├── cusum_depth_del.pl
│   ├── run_sample_base_depth_ratio.pl
│   ├── calc_depth_ratio.pl
│   ├── plot_depth_ratio.py
│   ├── extract_sa_split_reads.pl
│   ├── cluster_sa_split_reads.pl
│   ├── extract_valid_tlen.pl
│   ├── plot_tlen_distribution.py
│   ├── run_discordant_reads.pl
│   ├── merge_evidence.pl
│   ├── annotate_candidates.pl
│   ├── simulate_random_gene_wgs_del.pl
│   └── additional utility scripts
├── conf/
│   └── hcm_exondel.example.conf
├── db/
│   ├── exon_annotation_bed/
│   ├── canonical_transcripts/
│   ├── hcm_core_genes.txt
│   └── README
├── example/
├── test/
└── README.md
```

### Core scripts

| Script | Function |
|---|---|
| `HCMExonDel.pl` | Checks inputs and generates one sample-specific `run.sh`. |
| `extract_target_bam.pl` | Builds core-gene target regions and extracts a sorted, indexed target BAM. |
| `run_gene_mean_depth.pl` | Calculates gene depth, window depth, depth ratios, and consecutive low-depth windows. |
| `run_sample_cusum_depth.pl` | Runs CUSUM analysis for all retained per-gene depth files and combines the results. |
| `cusum_depth_del.pl` | Detects deletion-supporting intervals from base-level depth by CUSUM. |
| `run_sample_base_depth_ratio.pl` | Runs optional per-base depth-ratio calculations for all genes. |
| `calc_depth_ratio.pl` | Calculates per-base depth-ratio diagnostic output. |
| `plot_depth_ratio.py` | Generates a multi-page per-gene window-depth-ratio PDF. |
| `extract_sa_split_reads.pl` | Extracts breakpoint-supporting split/supplementary-alignment records. |
| `cluster_sa_split_reads.pl` | Clusters nearby split-read signals and filters them by unique-read support. |
| `extract_valid_tlen.pl` | Samples valid paired-end TLEN values for library QC. |
| `plot_tlen_distribution.py` | Plots the TLEN distribution and reports robust summary statistics. |
| `run_discordant_reads.pl` | Detects and clusters deletion-like discordant read pairs. |
| `merge_evidence.pl` | Uses split clusters as event backbones and evaluates depth/discordant support. |
| `annotate_candidates.pl` | Annotates candidate intervals against MANE RefSeq exons. |

---

## Input Files

### 1. Coordinate-sorted WGS BAM

Each sample must have a coordinate-sorted BAM file and a BAM index.

Accepted index names are:

```text
sample.bam.bai
sample.bai
```

Additional requirements:

- the BAM path in the BAM list must be absolute;
- the filename must end with `.bam`;
- the BAM reference build and chromosome naming must match the annotation files;
- sample names may contain letters, numbers, dots, underscores, and hyphens.

### 2. BAM list

The BAM list must contain exactly two TAB-delimited columns:

```text
SampleID<TAB>/absolute/path/to/sample.bam
```

Example:

```text
25B09089386 /data/project/HCM/25B09089386.final.merge.bam
25B09089387 /data/project/HCM/25B09089387.final.merge.bam
```

Rules:

- no header is required;
- blank lines are ignored;
- lines beginning with `#` are ignored;
- duplicate sample IDs are not allowed;
- duplicate BAM paths are not allowed.

### 3. MANE RefSeq exon annotation

Configured by:

```text
REFSEQ_MANE_SELECT_EXON_TXT
```

Required columns:

```text
Gene  Transcript  Exon  Chrom  Start  End  Strand
```

Coordinates must be **1-based closed**.

### 4. MANE RefSeq gene annotation

Configured by:

```text
REFSEQ_MANE_SELECT_GENE_TXT
```

Required columns:

```text
Gene  Transcript  Chrom  Start  End  Strand
```

Coordinates must be **1-based closed**.

The file is used to build target regions and perform depth analysis.

### 5. HCM core-gene list

Configured by:

```text
HCM_CORE_GENE_LIST
```

Example:

```text
Gene  Classification
MYBPC3  Definitive
MYH7  Definitive
FHOD3 Moderate
```

Only the first column is used as the gene symbol. A header named `Gene` is allowed. The classification column is optional.

A core gene absent from `REFSEQ_MANE_SELECT_GENE_TXT` is reported and skipped during target-region generation; it does not terminate the entire sample analysis.

---

## Configuration

Copy and edit the example configuration:

```bash
cp conf/hcm_exondel.example.conf conf/hcm_exondel.conf
```

The executable paths in the example file are environment-specific. In particular, update `PYTHON` and `SAMTOOLS` before running the pipeline on another system.

### Runtime and annotation

| Parameter | Example default | Description |
|---|---:|---|
| `PERL` | `/usr/bin/perl` | Perl executable. |
| `PYTHON` | environment-specific path | Python executable used for plotting. |
| `SAMTOOLS` | environment-specific path | samtools executable. |
| `REFSEQ_MANE_SELECT_EXON_TXT` | `db/exon_annotation_bed/RefSeq_MANE_Select.exon.txt` | MANE exon annotation. |
| `REFSEQ_MANE_SELECT_GENE_TXT` | `db/exon_annotation_bed/RefSeq_MANE_Select.gene.txt` | MANE gene annotation. |
| `HCM_CORE_GENE_LIST` | `db/hcm_core_genes.txt` | HCM core-gene list. |

Relative annotation paths are resolved against the project root.

### Target BAM and common read filters

| Parameter | Default | Description |
|---|---:|---|
| `TARGET_REGION_FLANK` | `5000` | Bases added on both sides of each core-gene interval when generating the target BAM. |
| `TARGET_BAM_THREADS` | `4` | Threads used during target BAM extraction/sorting. |
| `MIN_MAPQ` | `20` | Minimum mapping quality. |
| `EXCLUDE_DUPLICATES` | `1` | Exclude duplicate-marked reads when supported by the corresponding module. |

### Window-depth analysis

| Parameter | Default | Description |
|---|---:|---|
| `WINDOW_SIZE` | `500` | Window size in bp. |
| `WINDOW_STEP` | `200` | Step between adjacent windows in bp. |
| `MIN_WINDOW_SIZE` | `200` | Minimum retained window length in bp. |
| `MIN_GENE_MEAN_DEPTH` | `20` | Minimum gene mean depth for reliable ratio calculation. |
| `DEL_DEPTH_RATIO_CUTOFF` | `0.65` | Window depth-ratio threshold for deletion-supporting windows. |
| `MIN_CONSECUTIVE_DEL_WINDOWS` | `3` | Minimum number of consecutive deletion-supporting windows. |
| `KEEP_GENE_DEPTH_FILE` | `1` | Must remain `1`; CUSUM analysis depends on per-gene depth files. |
| `OUTPUT_BASE_DEPTH_RATIO` | `1` | Generate optional base-level ratio diagnostics. |
| `OUTPUT_BASE_DEPTH_RATIO_ALL` | `0` | Combine optional base-ratio outputs into sample-level files. Requires `OUTPUT_BASE_DEPTH_RATIO=1`. |

### CUSUM depth evidence

| Parameter | Default | Description |
|---|---:|---|
| `CUSUM_BASELINE` | `auto` | Baseline mode for the CUSUM calculation. |
| `CUSUM_BIN_SIZE` | `1` | CUSUM bin size in bp. |
| `CUSUM_K` | `0.1` | CUSUM reference/drift parameter. |
| `CUSUM_H` | `5` | CUSUM decision threshold. |
| `CUSUM_DEL_RATIO` | `0.65` | Ratio used to define deletion-like depth. |
| `CUSUM_MIN_BINS` | `3` | Minimum number of bins in a CUSUM interval. |
| `CUSUM_MIN_LEN` | `1` | Minimum interval length in bp. |
| `CUSUM_EDGE_RATIO` | `0.80` | Edge-extension ratio threshold. |
| `CUSUM_RECOVER_RATIO` | `0.80` | Recovery ratio threshold. |
| `CUSUM_RECOVER_BINS` | `3` | Consecutive recovery bins required to close an event. |

### Split-read clustering

| Parameter | Default | Description |
|---|---:|---|
| `SPLIT_READ_THREADS` | `4` | Threads used during split-read extraction. |
| `SA_SPLIT_CLUSTER_WINDOW` | `20` | Maximum breakpoint-coordinate distance for clustering. |
| `SA_SPLIT_MIN_SUPPORT_READS` | `5` | Minimum number of unique read names required for a passed cluster. |

### TLEN quality control

| Parameter | Default | Description |
|---|---:|---|
| `RUN_TLEN_QC` | `1` | Enable insert-size/TLEN QC. |
| `TLEN_MAX_PAIRS` | `1000000` | Maximum number of valid read pairs sampled. |
| `TLEN_THREADS` | `4` | Threads used for TLEN extraction. |
| `TLEN_RANDOM_SEED` | `20260701` | Random seed for reproducible sampling. |

The TLEN plot reports the median, MAD, robust standard deviation, and a diagnostic threshold of:

```text
median + 3 × 1.4826 × MAD
```

This diagnostic value does not automatically replace `MIN_DISCORDANT_INSERT_SIZE`. Tune the discordant threshold explicitly when library insert-size characteristics differ from the default setting.

### Discordant read-pair analysis

| Parameter | Default | Description |
|---|---:|---|
| `MIN_DISCORDANT_INSERT_SIZE` | `500` | Minimum absolute insert size for a discordant pair. |
| `MIN_DISCORDANT_READS` | `3` | Minimum supporting read pairs in a passed cluster. |
| `DISCORDANT_CLUSTER_DISTANCE` | `500` | Maximum distance used to group nearby discordant pairs. |
| `FILTER_DELETION_ORIENTATION` | `1` | Require deletion-like paired-end orientation. |
| `DISCORDANT_READ_THREADS` | `4` | Threads used during discordant-read extraction. |

### Evidence integration and annotation

| Parameter | Default | Description |
|---|---:|---|
| `EVIDENCE_OVERLAP_FRACTION` | `0.90` | Minimum fraction of the split interval that must be supported by qualifying depth or discordant evidence. |
| `MIN_EXON_OVERLAP_FRACTION` | `0.20` | Minimum fraction of an exon overlapped by a candidate for that exon to be reported as affected. |

### Runtime behavior

| Parameter | Default | Description |
|---|---:|---|
| `KEEP_TMP` | `0` | Remove the sample `tmp/` directory after successful completion. |
| `VERBOSE` | `1` | Enable verbose module output where supported. |

### Removed parameters

Do not include the following parameters in a current configuration file:

```text
ANALYZE_CORE_GENES_ONLY
CANDIDATE_BAM_FLANK
```

The main script explicitly rejects them.

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/Xiaohuaniu0032/HCMExonDel.git
cd HCMExonDel
```

### 2. Prepare the configuration

```bash
cp conf/hcm_exondel.example.conf conf/hcm_exondel.conf
vim conf/hcm_exondel.conf
```

At minimum, verify:

```text
PERL
PYTHON
SAMTOOLS
REFSEQ_MANE_SELECT_EXON_TXT
REFSEQ_MANE_SELECT_GENE_TXT
HCM_CORE_GENE_LIST
```

### 3. Prepare the BAM list

```bash
printf '25B09089386\t/absolute/path/to/25B09089386.final.merge.bam\n' > input_bam.list
```

### 4. Generate sample-specific shell scripts

```bash
perl HCMExonDel.pl \
  --config conf/hcm_exondel.conf \
  --bam-list input_bam.list \
  --outdir results
```

Expected output:

```text
results/25B09089386/25B09089386.run.sh
```

### 5. Run the analysis

```bash
bash results/25B09089386/25B09089386.run.sh
```

For multiple samples, execute or submit each generated `SAMPLE.run.sh` independently.

### 6. Existing output directory

By default, the main script stops when a sample output directory already exists. To allow writing into an existing sample directory:

```bash
perl HCMExonDel.pl \
  --config conf/hcm_exondel.conf \
  --bam-list input_bam.list \
  --outdir results \
  --force
```

For reproducible reruns, review or remove stale files in the existing sample directory before using `--force`.

---

## Output Structure

```text
OUTDIR/
└── SAMPLE/
    ├── SAMPLE.run.sh
    ├── 00.log/
    ├── 00.target_bam/
    ├── 01.depth/
    ├── 02.split_reads/
    ├── 03.discordant_reads/
    ├── 04.candidates/
    ├── 05.report/
    └── tmp/                      # removed when KEEP_TMP=0
```

### Log files

The generated workflow writes one log per step:

```text
00.extract_target_bam.log
01.gene_mean_depth.log
02.cusum_depth_evidence.log
03.base_depth_ratio.log
04.plot_depth_ratio.log
05.extract_sa_split_reads.log
06.cluster_sa_split_reads.log
07.extract_valid_tlen.log
08.plot_tlen_distribution.log
09.discordant_reads.log
10.merge_evidence.log
11.annotate_candidates.log
```

Skipped optional steps still receive a log containing the skip message.

---

## Key Output Files

### `00.target_bam/`

| File | Description |
|---|---|
| `SAMPLE.target.bam` | Sorted target BAM containing core genes plus configured flanks. |
| `SAMPLE.target.bam.bai` | BAM index. |
| `SAMPLE.target_regions.tsv` | Unmerged target-region records generated from the core-gene list and gene annotation. |
| `SAMPLE.target_regions.bed` | Merged BED intervals used for BAM extraction. |
| `SAMPLE.target_bam.summary.tsv` | Counts of requested, found, and skipped genes plus output paths and extraction settings. |

The BED file is generated internally; users do not need to provide a BED file.

### `01.depth/`

| File | Description | Used in final merge? |
|---|---|---:|
| `SAMPLE.depth_candidates.tsv` | Consecutive window-level low-depth candidate regions. | No |
| `SAMPLE.depth_candidates.all_window_ratio.tsv` | All analyzed windows and their depth ratios/status. | No |
| `SAMPLE.depth_candidates.del_windows.tsv` | Individual deletion-supporting windows. | No |
| `SAMPLE.depth_candidates.gene_mean_depth.tsv` | Gene-level mean-depth summary. | No |
| `SAMPLE.depth_candidates.window_depth.tsv` | Window-level depth summary. | No |
| `SAMPLE.depth_candidates.gene_depth_files/` | Per-gene base-depth files required by CUSUM. | Indirectly |
| `SAMPLE.depth_candidates.cusum/` | Per-gene CUSUM output. | Indirectly |
| `SAMPLE.depth_candidates.cusum.all.tsv` | Combined CUSUM deletion intervals. | **Yes** |
| `SAMPLE.depth_candidates.base_depth_ratio/` | Optional per-gene base-ratio diagnostics. | No |
| `SAMPLE.depth_candidates.base_depth_ratio.all.tsv` | Optional combined base-ratio output when enabled. | No |
| `SAMPLE.depth_candidates.base_depth.summary.all.tsv` | Optional combined base-ratio summary when enabled. | No |
| `SAMPLE.window_del.per_gene.pdf` | Multi-page core-gene window-depth-ratio plot. | No |

### `02.split_reads/`

| File | Description |
|---|---|
| `SAMPLE.split_reads.tsv` | Raw breakpoint-supporting split/supplementary-alignment records. |
| `SAMPLE.split_reads.clusters.tsv` | Passed split-read clusters; candidate backbones used by `merge_evidence.pl`. |

Additional diagnostic files may be produced by the clustering module depending on the current script implementation and filtering results.

### `03.discordant_reads/`

| File | Description |
|---|---|
| `SAMPLE.valid_tlen.tsv` | Sampled valid positive TLEN values when `RUN_TLEN_QC=1`. |
| `SAMPLE.tlen_distribution.pdf` | TLEN distribution plot. |
| `SAMPLE.tlen_distribution.pdf.summary.tsv` | TLEN summary statistics and robust diagnostic threshold. |
| `SAMPLE.discordant_reads.tsv` | Passed discordant-read clusters used by `merge_evidence.pl`. |

Additional raw/supporting/discarded-cluster files may be written by the discordant-read module.

### `04.candidates/`

| File | Description |
|---|---|
| `SAMPLE.merged_candidates.tsv` | All passed split-centered candidates with depth and discordant support metrics and an evidence level. |

There is no separate `SAMPLE.merged_candidates.tsv.all.tsv` in the current workflow.

### `05.report/`

| File | Description |
|---|---|
| `SAMPLE.annotated_candidates.tsv` | Final candidate report with MANE RefSeq exon annotation. |

---

## Candidate Interpretation

### Coordinates

Candidate coordinates are reported as **1-based closed intervals**. Candidate size is calculated as:

```text
End - Start + 1
```

The current merge step uses the split-read cluster to define:

```text
Best_Start
Best_End
Best_Size
Best_Evidence = Split
```

### Important merged-candidate columns

| Column | Description |
|---|---|
| `Cluster_ID` | Split-read cluster identifier. |
| `Split_Start`, `Split_End`, `Split_Size` | Candidate interval defined by the split cluster. |
| `Split_Reads` | Number of unique split-read names supporting the cluster. |
| `Split_Records` | Number of split-read records in the cluster. |
| `Depth_Support` | Whether CUSUM depth coverage reaches the overlap threshold. |
| `Depth_Covered_Bases` | Number of split-interval bases covered by qualifying depth intervals. |
| `Depth_Coverage_Fraction` | Covered fraction of the split interval. |
| `Depth_Record_Count` | Number of overlapping depth records contributing to support. |
| `Discordant_Read_Support` | Whether a qualifying discordant cluster supports the split interval. |
| `Discordant_Overlap_Fraction` | Best qualifying discordant overlap fraction. |
| `Discordant_Cluster_IDs` | Supporting discordant cluster identifiers. |
| `Discordant_Reads` | Number of supporting discordant read pairs. |
| `Median_Insert_Size` | Median insert size of the supporting discordant cluster. |
| `Evidence_Count` | Number of evidence types recorded for the candidate. |
| `Evidence_Level` | `High`, `Moderate`, or `Low`. |
| `Evidence_Types` | Evidence types supporting the candidate. |
| `Best_Start`, `Best_End`, `Best_Size` | Recommended split-defined candidate interval. |
| `Candidate_Status` | Candidate status assigned by the merge script. |
| `Comment` | Explanation of the evidence assignment. |

The pre-annotation `Gene` field may be `NA` because candidate intervals originate from split-read clusters. Use the final report's `Annotated_Gene`, `Annotated_Transcript`, and exon-overlap columns for gene/exon interpretation.

### Important annotation columns

| Column | Description |
|---|---|
| `Candidate_Region` | Final chromosome and candidate coordinates. |
| `Coordinate_Source` | Source used to select the final coordinates. |
| `Annotated_Gene` | Gene whose exon annotation passes the overlap requirement. |
| `Annotated_Transcript` | MANE transcript associated with the affected exon(s). |
| `Affected_Exons` | Exons meeting `MIN_EXON_OVERLAP_FRACTION`. |
| `Overlap_Exon_Count` | Number of affected exons. |
| `Fully_Covered_Exons` | Exons completely contained in the candidate interval. |
| `Partially_Overlapped_Exons` | Exons with qualifying partial overlap. |
| `Annotation_Status` | Whether a qualifying exon overlap was found. |
| `Exon_Overlap_Detail` | Exon coordinates, overlap length, and overlap fraction. |

A `High` evidence level describes sequencing evidence for a split-centered event. It does **not** by itself mean that the event affects a coding exon. Always interpret `Evidence_Level` together with `Annotation_Status` and the exon-overlap columns.

---

## Test Example: 25B09089386

The repository contains an example result directory:

```text
test/test_results/25B09089386/
```

### Target-region extraction

According to `00.target_bam/25B09089386.target_bam.summary.tsv`:

| Metric | Result |
|---|---:|
| Core genes in list | 35 |
| Core genes found in MANE gene annotation | 34 |
| Core genes skipped | 1 |
| Skipped gene | `MT-TI` |
| Raw gene regions | 34 |
| Merged target BED regions | 34 |
| Target-region flank | 5,000 bp |

The skipped gene is reported in the summary rather than causing the sample workflow to terminate.

### Merged candidates

`04.candidates/25B09089386.merged_candidates.tsv` contains seven split-centered candidates:

| Evidence level | Count |
|---|---:|
| `High` | 3 |
| `Moderate` | 2 |
| `Low` | 2 |
| **Total** | **7** |

### Exon annotation

Only one of the seven candidate intervals meets the configured MANE exon-overlap requirement:

| Field | Result |
|---|---|
| Sample | `25B09089386` |
| Gene | `FHOD3` |
| Transcript | `NM_001281740` |
| Region | `chr18:34232240-34241309` |
| Size | 9,070 bp |
| Split-read support | 16 unique reads |
| Depth support | Yes |
| Depth coverage fraction | 0.9994 |
| Discordant support | Yes |
| Evidence level | `High` |
| Affected exons | `EX12`, `EX13`, `EX14` |
| Fully covered exons | `EX12`, `EX13`, `EX14` |

The other six intervals are retained as sequencing-evidence candidates but are reported as having no qualifying exon overlap. This distinction is intentional:

```text
Evidence strength != exon consequence
```

The test dataset demonstrates expected file generation and candidate prioritization. It should not be interpreted as a formal sensitivity, specificity, or clinical-performance validation dataset.

---

## Manual Review and Validation

All candidates should be reviewed before biological or clinical interpretation.

Recommended review procedure:

1. Open the original WGS BAM and/or `SAMPLE.target.bam` in IGV.
2. Inspect depth across the candidate interval and adjacent exons.
3. Confirm that the depth decrease is consistent across the proposed deletion.
4. Inspect split reads and supplementary alignments near both breakpoints.
5. Inspect discordant read pairs spanning the interval.
6. Check mapping quality, duplicate status, repetitive sequence, segmental duplications, and nearby alignment artifacts.
7. Confirm that chromosome names and reference coordinates match the annotation build.
8. Compare affected and unaffected family members when available.
9. Prioritize candidates affecting coding exons and segregating with the phenotype.

Possible orthogonal validation methods include:

- breakpoint PCR;
- Sanger sequencing across the breakpoint;
- ddPCR;
- MLPA;
- qPCR;
- RNA-level confirmation when relevant and biologically interpretable;
- family segregation analysis.

---

## Running Individual Modules

The recommended approach is to generate and execute `SAMPLE.run.sh`. The following commands illustrate the current module interfaces. Downstream modules should normally use `SAMPLE.target.bam`, not the original WGS BAM.

### Target BAM extraction

```bash
perl bin/extract_target_bam.pl \
  --config conf/hcm_exondel.conf \
  --sample SAMPLE \
  --bam /absolute/path/to/SAMPLE.bam \
  --outdir results/SAMPLE/00.target_bam
```

### Window-depth analysis

```bash
perl bin/run_gene_mean_depth.pl \
  --config conf/hcm_exondel.conf \
  --bam results/SAMPLE/00.target_bam/SAMPLE.target.bam \
  --sample SAMPLE \
  --out results/SAMPLE/01.depth/SAMPLE.depth_candidates.tsv
```

### Sample-level CUSUM analysis

```bash
perl bin/run_sample_cusum_depth.pl \
  --config conf/hcm_exondel.conf \
  --in-dir results/SAMPLE/01.depth/SAMPLE.depth_candidates.gene_depth_files \
  --out-dir results/SAMPLE/01.depth/SAMPLE.depth_candidates.cusum \
  --out results/SAMPLE/01.depth/SAMPLE.depth_candidates.cusum.all.tsv \
  --perl /usr/bin/perl \
  --script bin/cusum_depth_del.pl
```

### Depth-ratio plotting

```bash
python3 bin/plot_depth_ratio.py \
  --input results/SAMPLE/01.depth/SAMPLE.depth_candidates.all_window_ratio.tsv \
  --output results/SAMPLE/01.depth/SAMPLE.window_del.per_gene.pdf
```

### Split-read extraction

```bash
perl bin/extract_sa_split_reads.pl \
  --conf conf/hcm_exondel.conf \
  --bam results/SAMPLE/00.target_bam/SAMPLE.target.bam \
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
  --bam results/SAMPLE/00.target_bam/SAMPLE.target.bam \
  --sample SAMPLE \
  --out results/SAMPLE/03.discordant_reads/SAMPLE.discordant_reads.tsv
```

### Evidence merging

```bash
perl bin/merge_evidence.pl \
  --sample SAMPLE \
  --depth results/SAMPLE/01.depth/SAMPLE.depth_candidates.cusum.all.tsv \
  --split results/SAMPLE/02.split_reads/SAMPLE.split_reads.clusters.tsv \
  --discordant results/SAMPLE/03.discordant_reads/SAMPLE.discordant_reads.tsv \
  --min-overlap 0.90 \
  --out results/SAMPLE/04.candidates/SAMPLE.merged_candidates.tsv
```

Important points:

- `merge_evidence.pl` does not read the configuration file directly;
- pass `SAMPLE.depth_candidates.cusum.all.tsv`, not the window-level `SAMPLE.depth_candidates.tsv`;
- pass the clustered split-read file, not the raw split-read file;
- the main script reads `EVIDENCE_OVERLAP_FRACTION` and supplies it as `--min-overlap`.

### Candidate annotation

```bash
perl bin/annotate_candidates.pl \
  --config conf/hcm_exondel.conf \
  --input results/SAMPLE/04.candidates/SAMPLE.merged_candidates.tsv \
  --out results/SAMPLE/05.report/SAMPLE.annotated_candidates.tsv
```

---

## Troubleshooting

### BAM index not found

```text
[ERROR] BAM index not found
```

Create an index:

```bash
samtools index sample.bam
```

### BAM index is older than the BAM

A stale index can cause incorrect regional extraction or warnings. Rebuild it:

```bash
rm -f sample.bam.bai sample.bai
samtools index sample.bam
```

### Invalid BAM list format

The BAM list must contain exactly two TAB-delimited columns. A sequence of spaces is not equivalent to a TAB.

Check hidden characters with:

```bash
cat -A input_bam.list
```

### Output directory already exists

Use `--force` only after confirming that writing into the existing directory is intended. Remove stale intermediate files when parameter settings or code versions have changed.

### Deprecated parameter error

Remove either of these obsolete settings from the configuration:

```text
ANALYZE_CORE_GENES_ONLY
CANDIDATE_BAM_FLANK
```

### Core gene not found in the gene annotation

The target-BAM summary lists skipped genes. Check:

- the gene symbol spelling;
- whether the symbol is present in the MANE gene annotation;
- whether mitochondrial or non-MANE genes require a separate annotation strategy;
- whether the gene list contains aliases rather than the annotation's approved symbol.

A missing gene is skipped and reported; it does not invalidate successfully resolved genes.

### Chromosome naming mismatch

Examples of incompatible naming include:

```text
chr18 vs 18
chrM  vs MT
```

Ensure that the BAM and annotation files use the same reference build and chromosome naming convention.

### `KEEP_GENE_DEPTH_FILE=0` error

The current CUSUM step requires per-gene base-depth files. Keep:

```text
KEEP_GENE_DEPTH_FILE=1
```

### `OUTPUT_BASE_DEPTH_RATIO_ALL=1` error

Combined base-ratio output requires base-ratio generation:

```text
OUTPUT_BASE_DEPTH_RATIO=1
OUTPUT_BASE_DEPTH_RATIO_ALL=1
```

### No split-read candidates

Because the current merge is split-read-centered, no passed split cluster means no merged candidate, even when depth or discordant signals exist.

Possible causes include:

- a true breakpoint in repetitive or poorly mappable sequence;
- missing or filtered supplementary alignments;
- insufficient breakpoint-spanning reads;
- `SA_SPLIT_MIN_SUPPORT_READS` being too strict;
- low mapping quality;
- breakpoint locations falling outside the extracted target region and flank.

Review depth and discordant diagnostic outputs separately when investigating such cases.

### Candidate remains `Low` despite depth overlap

Depth support requires coverage of at least `EVIDENCE_OVERLAP_FRACTION` of the split interval. Partial overlap below the threshold is recorded but does not qualify as depth support.

Review:

```text
Depth_Covered_Bases
Depth_Coverage_Fraction
Depth_Range
```

### Weak discordant support

Possible causes include:

- the deletion is smaller than `MIN_DISCORDANT_INSERT_SIZE`;
- the library insert-size distribution differs from the configured threshold;
- too few supporting read pairs;
- orientation filtering;
- mapping ambiguity around the breakpoints.

Use the TLEN summary as a diagnostic guide, then tune `MIN_DISCORDANT_INSERT_SIZE` explicitly when justified.

### Empty merged output

An empty merged file is a valid outcome when no split cluster passes filtering. The script should retain the header so downstream annotation can terminate cleanly.

Inspect:

```text
02.split_reads/SAMPLE.split_reads.tsv
02.split_reads/SAMPLE.split_reads.clusters.tsv
00.log/05.extract_sa_split_reads.log
00.log/06.cluster_sa_split_reads.log
```

### High-evidence candidate has no exon annotation

`High` means that split, depth, and discordant evidence support the genomic event. It does not guarantee exon overlap.

Inspect:

```text
Annotation_Status
Annotated_Gene
Affected_Exons
Exon_Overlap_Detail
```

Such a candidate may be intronic, intergenic, below the exon-overlap threshold, affected by annotation/reference mismatch, or artifactual.

---

## Limitations

- HCMExonDel is designed for WGS BAM data and is not optimized for WES or targeted-panel data.
- The pipeline analyzes only genes listed in `HCM_CORE_GENE_LIST`.
- The final merge is split-read-centered; events without a passed split-read cluster are not reported as merged candidates.
- Breakpoints in repetitive, low-complexity, or poorly mappable regions may lack reliable split-read support.
- Read-depth analysis is affected by sequencing depth, GC bias, mappability, alignment quality, and local coverage variability.
- Small deletions may not produce sufficient depth or discordant-pair evidence.
- Large deletions may extend beyond the target gene plus configured flank, reducing breakpoint evidence in the target BAM.
- A fixed `MIN_DISCORDANT_INSERT_SIZE` may not be optimal for every sequencing library.
- Exon annotation depends on the selected MANE transcript file and `MIN_EXON_OVERLAP_FRACTION`.
- Evidence level is not equivalent to pathogenicity, clinical classification, or exon consequence.
- The included test sample demonstrates workflow behavior but does not establish analytical sensitivity, specificity, precision, or limit of detection.
- All candidates require manual review and independent experimental validation.

---

## Citation

A formal publication citation is not yet available. When using HCMExonDel in a report or manuscript, cite the repository and describe the analysis version and configuration parameters.

Suggested methods text:

```text
Candidate exon-level deletions in hypertrophic cardiomyopathy-related core genes were detected from whole-genome sequencing BAM files using HCMExonDel. The pipeline generated a target BAM spanning predefined core genes and flanking regions, identified split-read clusters as candidate event backbones, evaluated supporting base-level CUSUM depth and discordant read-pair evidence, and annotated candidate intervals against MANE RefSeq exon coordinates. All candidates were manually reviewed and required independent validation.
```

Also report, at minimum:

- repository commit or release;
- reference genome build;
- MANE annotation version;
- HCM core-gene list version;
- principal depth, split-read, discordant-read, and overlap thresholds.

---

## License

No license file is currently included in the repository. Add an explicit license before distributing or reusing the software beyond the intended project context.

---

## Contact

Maintainer: `Xiaohuaniu0032`

Repository:

```text
https://github.com/Xiaohuaniu0032/HCMExonDel
```

---

## Disclaimer

HCMExonDel is provided for research use only. It is not intended for direct clinical diagnosis. Candidate deletion events must be reviewed by qualified personnel and validated using independent methods before clinical interpretation or reporting.


