#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot exon/CDS/UTR structures for HCM core genes based on MANE Select transcripts
and hg19 ncbiRefSeq GTF annotation.

Input:
  1. hcm_core_genes.txt
     Required columns:
       Gene    Classification

  2. RefSeq_MANE_Select.xls
     Required columns:
       name    RefSeq_prot    MANE_status
     Note:
       Although the file suffix is .xls, it may be a tab-delimited text file.

  3. hg19.ncbiRefSeq.gtf
     Required attributes:
       gene_id, transcript_id, gene_name, exon_number

Output:
  One PDF file, one gene per page.
"""

import argparse
import os
import re
import sys

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter


def strip_version(transcript_id):
    """
    Remove transcript version.

    Example:
      NM_001281739.3 -> NM_001281739
    """
    if pd.isna(transcript_id):
        return transcript_id

    return str(transcript_id).split(".")[0]


def read_table_auto(path):
    """
    Read tab-delimited text or real Excel file automatically.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")

    # Most .xls files in this workflow are tab-delimited text files.
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)

        if df.shape[1] >= 2:
            return df
    except Exception:
        pass

    # Fallback for real Excel files.
    try:
        return pd.read_excel(path, dtype=str)
    except Exception as e:
        raise RuntimeError(f"Failed to read file as TSV or Excel: {path}\n{e}")


def read_core_genes(core_gene_file):
    """
    Read HCM core gene list and deduplicate genes while preserving order.
    """
    df = read_table_auto(core_gene_file)

    if "Gene" not in df.columns:
        raise ValueError(f"'Gene' column not found in {core_gene_file}")

    genes = (
        df["Gene"]
        .dropna()
        .astype(str)
        .str.strip()
    )

    genes = [g for g in genes if g and g.lower() != "gene"]

    seen = set()
    unique_genes = []

    for gene in genes:
        if gene not in seen:
            unique_genes.append(gene)
            seen.add(gene)

    return unique_genes


def read_mane_transcripts(mane_file):
    """
    Read MANE Select transcript table.

    Expected columns:
      name
      RefSeq_prot
      MANE_status
    """
    df = read_table_auto(mane_file)

    required_cols = {"name", "RefSeq_prot"}
    missing = required_cols - set(df.columns)

    if missing:
        raise ValueError(
            f"Missing required columns in {mane_file}: {', '.join(sorted(missing))}"
        )

    if "MANE_status" in df.columns:
        df = df[df["MANE_status"].astype(str).str.contains("MANE Select", na=False)]

    mane = {}

    for _, row in df.iterrows():
        gene = str(row["name"]).strip()
        tx = str(row["RefSeq_prot"]).strip()

        if not gene or gene.lower() == "nan":
            continue

        if not tx or tx.lower() == "nan":
            continue

        mane[gene] = strip_version(tx)

    return mane


def parse_gtf_attributes(attr_string):
    """
    Parse GTF attributes.

    Example:
      gene_id "FHOD3"; transcript_id "NM_001281739.3"; gene_name "FHOD3";
    """
    attrs = {}

    for item in attr_string.strip().split(";"):
        item = item.strip()

        if not item:
            continue

        m = re.match(r'(\S+)\s+"([^"]+)"', item)

        if m:
            attrs[m.group(1)] = m.group(2)

    return attrs


def parse_gtf_for_targets(gtf_file, target_gene_to_tx):
    """
    Parse GTF and keep records matching target genes and MANE transcripts.
    """
    gene_records = {}

    target_genes = set(target_gene_to_tx.keys())

    with open(gtf_file, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs_raw = fields

            try:
                start = int(start)
                end = int(end)
            except ValueError:
                continue

            attrs = parse_gtf_attributes(attrs_raw)

            gene_name = attrs.get("gene_name") or attrs.get("gene_id")
            transcript_id_full = attrs.get("transcript_id")

            if not gene_name or not transcript_id_full:
                continue

            transcript_id = strip_version(transcript_id_full)

            if gene_name not in target_genes:
                continue

            expected_tx = target_gene_to_tx[gene_name]

            if transcript_id != expected_tx:
                continue

            if gene_name not in gene_records:
                gene_records[gene_name] = {
                    "transcript": transcript_id,
                    "transcript_full": transcript_id_full,
                    "chrom": chrom,
                    "strand": strand,
                    "tx_start": None,
                    "tx_end": None,
                    "exons": [],
                    "cds": [],
                    "utr5": [],
                    "utr3": [],
                }

            rec = gene_records[gene_name]

            if feature == "transcript":
                rec["tx_start"] = start
                rec["tx_end"] = end

            elif feature == "exon":
                exon_number = attrs.get("exon_number", "NA")
                exon_id = attrs.get("exon_id", f"{transcript_id_full}.{exon_number}")

                rec["exons"].append(
                    {
                        "start": start,
                        "end": end,
                        "exon_number": exon_number,
                        "exon_id": exon_id,
                    }
                )

            elif feature == "CDS":
                rec["cds"].append(
                    {
                        "start": start,
                        "end": end,
                    }
                )

            elif feature == "5UTR":
                rec["utr5"].append(
                    {
                        "start": start,
                        "end": end,
                    }
                )

            elif feature == "3UTR":
                rec["utr3"].append(
                    {
                        "start": start,
                        "end": end,
                    }
                )

    for gene, rec in gene_records.items():
        rec["exons"].sort(key=lambda x: x["start"])
        rec["cds"].sort(key=lambda x: x["start"])
        rec["utr5"].sort(key=lambda x: x["start"])
        rec["utr3"].sort(key=lambda x: x["start"])

        if rec["tx_start"] is None or rec["tx_end"] is None:
            if rec["exons"]:
                rec["tx_start"] = min(x["start"] for x in rec["exons"])
                rec["tx_end"] = max(x["end"] for x in rec["exons"])

    return gene_records


def infer_utr_from_exon_cds(exons, cds):
    """
    If GTF does not contain explicit 5UTR / 3UTR records,
    infer non-CDS parts within exons.

    This returns generic UTR intervals, not distinguishing 5UTR and 3UTR.
    """
    inferred = []

    if not exons or not cds:
        return inferred

    for exon in exons:
        exon_start = exon["start"]
        exon_end = exon["end"]

        covered = []

        for c in cds:
            s = max(exon_start, c["start"])
            e = min(exon_end, c["end"])

            if s <= e:
                covered.append((s, e))

        if not covered:
            inferred.append({"start": exon_start, "end": exon_end})
            continue

        covered.sort()

        cursor = exon_start

        for s, e in covered:
            if cursor < s:
                inferred.append({"start": cursor, "end": s - 1})

            cursor = max(cursor, e + 1)

        if cursor <= exon_end:
            inferred.append({"start": cursor, "end": exon_end})

    return inferred


def add_interval_patch(ax, start, end, y, height, color):
    """
    Draw one interval without black border.
    """
    width = end - start + 1

    rect = Rectangle(
        (start, y - height / 2),
        width,
        height,
        facecolor=color,
        edgecolor="none",
        linewidth=0,
    )

    ax.add_patch(rect)


def format_full_integer(x, pos):
    """
    Force genomic coordinate labels to be full integers.
    """
    return f"{int(x):d}"


def get_first_last_exon_numbers(exons):
    """
    Get the first and last exon numbers based on genomic order.

    Since exons have already been sorted by genomic start,
    this function labels the first and last exon shown in the figure.
    """
    if not exons:
        return set()

    first_exon_number = str(exons[0].get("exon_number", "NA"))
    last_exon_number = str(exons[-1].get("exon_number", "NA"))

    return {first_exon_number, last_exon_number}


def plot_gene_structure(ax, gene, rec):
    chrom = rec["chrom"]
    strand = rec["strand"]
    transcript = rec["transcript_full"]
    tx_start = rec["tx_start"]
    tx_end = rec["tx_end"]

    exons = rec["exons"]
    cds = rec["cds"]
    utr5 = rec["utr5"]
    utr3 = rec["utr3"]

    if not utr5 and not utr3:
        inferred_utr = infer_utr_from_exon_cds(exons, cds)
    else:
        inferred_utr = []

    y = 0

    # Transcript backbone
    ax.hlines(y, tx_start, tx_end, color="black", linewidth=0.8)

    # Exon background
    for exon in exons:
        add_interval_patch(
            ax,
            exon["start"],
            exon["end"],
            y,
            height=0.26,
            color="#D9D9D9",
        )

    # 5' UTR
    for u in utr5:
        add_interval_patch(
            ax,
            u["start"],
            u["end"],
            y,
            height=0.22,
            color="#66C2A5",
        )

    # 3' UTR
    for u in utr3:
        add_interval_patch(
            ax,
            u["start"],
            u["end"],
            y,
            height=0.22,
            color="#FC8D62",
        )

    # Inferred UTR if no explicit UTR feature is available
    for u in inferred_utr:
        add_interval_patch(
            ax,
            u["start"],
            u["end"],
            y,
            height=0.22,
            color="#B3B3B3",
        )

    # CDS
    for c in cds:
        add_interval_patch(
            ax,
            c["start"],
            c["end"],
            y,
            height=0.50,
            color="#4C72B0",
        )

    # Exon labels: only label the first and last exon number
    exon_numbers_to_label = get_first_last_exon_numbers(exons)

    for exon in exons:
        exon_number = str(exon.get("exon_number", "NA"))

        if exon_number not in exon_numbers_to_label:
            continue

        exon_mid = (exon["start"] + exon["end"]) / 2

        ax.text(
            exon_mid,
            y + 0.38,
            exon_number,
            ha="center",
            va="bottom",
            fontsize=7,
            rotation=90,
        )

    # Direction arrow
    arrow_y = -0.45

    if strand == "+":
        ax.annotate(
            "",
            xy=(tx_end, arrow_y),
            xytext=(tx_start, arrow_y),
            arrowprops=dict(arrowstyle="->", linewidth=0.8),
        )
    else:
        ax.annotate(
            "",
            xy=(tx_start, arrow_y),
            xytext=(tx_end, arrow_y),
            arrowprops=dict(arrowstyle="->", linewidth=0.8),
        )

    ax.text(
        tx_start,
        -0.75,
        f"{chrom}:{tx_start:,}-{tx_end:,} ({strand})",
        ha="left",
        va="top",
        fontsize=8,
    )

    # Title: gene | transcript
    ax.set_title(
        f"{gene} | {transcript}",
        fontsize=12,
        fontweight="bold",
    )

    span = tx_end - tx_start + 1
    pad = max(int(span * 0.03), 100)

    ax.set_xlim(tx_start - pad, tx_end + pad)
    ax.set_ylim(-1.2, 1.2)

    ax.set_yticks([])
    ax.set_xlabel("Genomic position, hg19", fontsize=9)

    # Full integer coordinate labels.
    # Do not use ax.ticklabel_format() after FuncFormatter,
    # because ticklabel_format() only works with ScalarFormatter.
    ax.xaxis.set_major_formatter(FuncFormatter(format_full_integer))
    ax.get_xaxis().get_offset_text().set_visible(False)

    for label in ax.get_xticklabels():
        label.set_rotation(30)
        label.set_ha("right")
        label.set_fontsize(7)

    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    # Legend without black border
    legend_items = [
        Rectangle((0, 0), 1, 1, facecolor="#D9D9D9", edgecolor="none", label="Exon"),
        Rectangle((0, 0), 1, 1, facecolor="#4C72B0", edgecolor="none", label="CDS"),
        Rectangle((0, 0), 1, 1, facecolor="#66C2A5", edgecolor="none", label="5' UTR"),
        Rectangle((0, 0), 1, 1, facecolor="#FC8D62", edgecolor="none", label="3' UTR"),
    ]

    if inferred_utr:
        legend_items.append(
            Rectangle(
                (0, 0),
                1,
                1,
                facecolor="#B3B3B3",
                edgecolor="none",
                label="Inferred UTR",
            )
        )

    ax.legend(
        handles=legend_items,
        loc="upper right",
        fontsize=7,
        frameon=False,
        ncol=4,
    )


def write_missing_report(path, missing_mane, missing_gtf):
    with open(path, "w") as out:
        out.write("Category\tGene\tReason\n")

        for gene in missing_mane:
            out.write(
                f"Missing_MANE\t{gene}\tGene not found in RefSeq_MANE_Select file\n"
            )

        for gene in missing_gtf:
            out.write(
                f"Missing_GTF\t{gene}\tMANE transcript not found in GTF\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description="Plot exon/CDS/UTR structures for HCM core genes using MANE Select transcripts."
    )

    parser.add_argument(
        "--core-genes",
        required=True,
        help="Input hcm_core_genes.txt file. Required column: Gene",
    )

    parser.add_argument(
        "--mane",
        required=True,
        help="Input RefSeq_MANE_Select.xls file. Required columns: name, RefSeq_prot",
    )

    parser.add_argument(
        "--gtf",
        required=True,
        help="Input hg19 ncbiRefSeq GTF file",
    )

    parser.add_argument(
        "--out-pdf",
        default="HCM_core_gene_structures.pdf",
        help="Output PDF file. Default: HCM_core_gene_structures.pdf",
    )

    parser.add_argument(
        "--missing-report",
        default="HCM_core_gene_missing_report.tsv",
        help="Output missing gene report. Default: HCM_core_gene_missing_report.tsv",
    )

    args = parser.parse_args()

    print("[INFO] Reading HCM core genes...", file=sys.stderr)
    core_genes = read_core_genes(args.core_genes)
    print(f"[INFO] Unique core genes: {len(core_genes)}", file=sys.stderr)

    print("[INFO] Reading MANE Select transcript table...", file=sys.stderr)
    mane = read_mane_transcripts(args.mane)
    print(f"[INFO] MANE records loaded: {len(mane)}", file=sys.stderr)

    target_gene_to_tx = {}
    missing_mane = []

    for gene in core_genes:
        if gene not in mane:
            missing_mane.append(gene)
            print(f"[WARN] MANE transcript not found for gene: {gene}", file=sys.stderr)
            continue

        target_gene_to_tx[gene] = mane[gene]

    print(f"[INFO] Genes with MANE transcript: {len(target_gene_to_tx)}", file=sys.stderr)

    print("[INFO] Parsing GTF annotation...", file=sys.stderr)
    gene_records = parse_gtf_for_targets(args.gtf, target_gene_to_tx)
    print(f"[INFO] Genes found in GTF: {len(gene_records)}", file=sys.stderr)

    missing_gtf = []

    for gene in target_gene_to_tx:
        if gene not in gene_records:
            missing_gtf.append(gene)
            print(
                f"[WARN] GTF record not found for gene/transcript: "
                f"{gene} / {target_gene_to_tx[gene]}",
                file=sys.stderr,
            )

    if not gene_records:
        raise RuntimeError(
            "No gene structure records found. Please check gene names, transcript IDs, and GTF file."
        )

    print(f"[INFO] Writing PDF: {args.out_pdf}", file=sys.stderr)

    with PdfPages(args.out_pdf) as pdf:
        for gene in core_genes:
            if gene not in gene_records:
                continue

            rec = gene_records[gene]

            if not rec["exons"]:
                print(f"[WARN] No exon records for gene: {gene}", file=sys.stderr)
                continue

            fig, ax = plt.subplots(figsize=(12, 3.2))
            plot_gene_structure(ax, gene, rec)

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    write_missing_report(args.missing_report, missing_mane, missing_gtf)

    print("[INFO] Done.", file=sys.stderr)
    print(f"[INFO] Output PDF: {args.out_pdf}", file=sys.stderr)
    print(f"[INFO] Missing report: {args.missing_report}", file=sys.stderr)


if __name__ == "__main__":
    main()

    