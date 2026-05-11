#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_depth_ratio.py

Plot window-based depth ratio for HCMExonDel.

Input:
    *.all_window_ratio.tsv from run_gene_mean_depth.pl

Required columns:
    SampleID
    Gene
    Transcript
    Region
    Window_ID
    Chrom
    Start
    End
    Length
    Region_Mean_Depth
    Window_Mean_Depth
    Depth_Ratio
    Window_Status
    Comment

Example:
    python scripts/plot_depth_ratio.py \
        --input results/SAMPLE001/01.depth/SAMPLE001.depth_candidates.all_window_ratio.tsv \
        --out results/SAMPLE001/05.report/SAMPLE001.depth_ratio.png
"""

import argparse
import os
import sys
import math
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot window-based depth ratio from HCMExonDel output."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input *.all_window_ratio.tsv file."
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Output figure path. Supported: png, pdf, svg, tiff."
    )

    parser.add_argument(
        "--gene",
        default=None,
        help="Plot only one gene. Default: plot all genes in separate figures if multiple genes exist."
    )

    parser.add_argument(
        "--chrom",
        default=None,
        help="Plot only one chromosome."
    )

    parser.add_argument(
        "--region",
        default=None,
        help="Plot only one region name."
    )

    parser.add_argument(
        "--ratio-cutoff",
        type=float,
        default=0.65,
        help="Depth ratio cutoff line. Default: 0.65."
    )

    parser.add_argument(
        "--het-low",
        type=float,
        default=0.35,
        help="Lower bound of heterozygous deletion-like ratio. Default: 0.35."
    )

    parser.add_argument(
        "--het-high",
        type=float,
        default=0.70,
        help="Upper bound of heterozygous deletion-like ratio. Default: 0.70."
    )

    parser.add_argument(
        "--fig-width",
        type=float,
        default=12,
        help="Figure width. Default: 12."
    )

    parser.add_argument(
        "--fig-height",
        type=float,
        default=4.5,
        help="Figure height. Default: 4.5."
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Figure DPI. Default: 300."
    )

    parser.add_argument(
        "--title",
        default=None,
        help="Custom figure title."
    )

    return parser.parse_args()


def check_required_columns(df, input_file):
    required = [
        "SampleID",
        "Gene",
        "Transcript",
        "Region",
        "Window_ID",
        "Chrom",
        "Start",
        "End",
        "Length",
        "Region_Mean_Depth",
        "Window_Mean_Depth",
        "Depth_Ratio",
        "Window_Status",
    ]

    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError(
            f"Required columns missing in {input_file}: {','.join(missing)}"
        )


def prepare_dataframe(df):
    df = df.copy()

    numeric_cols = [
        "Start",
        "End",
        "Length",
        "Region_Mean_Depth",
        "Window_Mean_Depth",
        "Depth_Ratio",
    ]

    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["Start", "End", "Depth_Ratio"])

    df["Midpoint"] = (df["Start"] + df["End"]) / 2

    df = df.sort_values(["Chrom", "Gene", "Region", "Start", "End"])

    return df


def is_low_depth_status(status):
    if pd.isna(status):
        return False

    status = str(status)

    low_status = {
        "Low_depth_window",
        "Very_low_depth_window",
        "Mild_low_depth_window",
    }

    return status in low_status


def get_low_depth_segments(df):
    """
    Merge consecutive low-depth windows for visual shading.
    """
    segments = []

    sub = df.sort_values(["Start", "End"]).copy()
    sub["Is_Low"] = sub["Window_Status"].apply(is_low_depth_status)

    current_start = None
    current_end = None

    for _, row in sub.iterrows():
        if row["Is_Low"]:
            if current_start is None:
                current_start = row["Start"]
                current_end = row["End"]
            else:
                if row["Start"] <= current_end + 1:
                    current_end = max(current_end, row["End"])
                else:
                    segments.append((current_start, current_end))
                    current_start = row["Start"]
                    current_end = row["End"]
        else:
            if current_start is not None:
                segments.append((current_start, current_end))
                current_start = None
                current_end = None

    if current_start is not None:
        segments.append((current_start, current_end))

    return segments


def format_bp(x):
    try:
        x = float(x)
    except Exception:
        return str(x)

    if abs(x) >= 1_000_000:
        return f"{x / 1_000_000:.2f} Mb"
    if abs(x) >= 1_000:
        return f"{x / 1_000:.1f} kb"

    return str(int(x))


def plot_one_group(df, args, output_path):
    sample = str(df["SampleID"].iloc[0])
    gene = str(df["Gene"].iloc[0])
    transcript = str(df["Transcript"].iloc[0])
    region = str(df["Region"].iloc[0])
    chrom = str(df["Chrom"].iloc[0])

    df = df.sort_values(["Start", "End"]).copy()

    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))

    ax.plot(
        df["Midpoint"],
        df["Depth_Ratio"],
        marker="o",
        linewidth=1,
        markersize=3,
        label="Window depth ratio"
    )

    low_df = df[df["Window_Status"].apply(is_low_depth_status)]
    if not low_df.empty:
        ax.scatter(
            low_df["Midpoint"],
            low_df["Depth_Ratio"],
            s=18,
            label="Low-depth windows"
        )

    ax.axhline(
        y=args.ratio_cutoff,
        linestyle="--",
        linewidth=1,
        label=f"Ratio cutoff = {args.ratio_cutoff}"
    )

    ax.axhline(
        y=args.het_low,
        linestyle=":",
        linewidth=1,
        label=f"Het-del lower = {args.het_low}"
    )

    ax.axhline(
        y=args.het_high,
        linestyle=":",
        linewidth=1,
        label=f"Het-del upper = {args.het_high}"
    )

    segments = get_low_depth_segments(df)
    for start, end in segments:
        ax.axvspan(start, end, alpha=0.15)

    if args.title:
        title = args.title
    else:
        title = f"{sample} | {gene} | {transcript} | {chrom}:{int(df['Start'].min())}-{int(df['End'].max())}"

    ax.set_title(title)
    ax.set_xlabel(f"Genomic position on {chrom}")
    ax.set_ylabel("Window depth ratio")

    ymax = max(1.5, df["Depth_Ratio"].max() * 1.2)
    if math.isfinite(ymax):
        ax.set_ylim(0, ymax)

    xmin = df["Start"].min()
    xmax = df["End"].max()
    ax.set_xlim(xmin, xmax)

    xticks = ax.get_xticks()
    ax.set_xticklabels([format_bp(x) for x in xticks], rotation=30, ha="right")

    subtitle = f"Region: {region}; low-depth segments: {len(segments)}"
    ax.text(
        0.01,
        0.98,
        subtitle,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9
    )

    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()

    ext = os.path.splitext(output_path)[1].lower()

    if ext in [".tif", ".tiff"]:
        fig.savefig(output_path, dpi=args.dpi, pil_kwargs={"compression": "tiff_lzw"})
    else:
        fig.savefig(output_path, dpi=args.dpi)

    plt.close(fig)


def build_output_path(base_out, gene, region, index, total):
    if total == 1:
        return base_out

    root, ext = os.path.splitext(base_out)
    safe_gene = str(gene).replace("/", "_").replace("|", "_").replace(":", "_")
    safe_region = str(region).replace("/", "_").replace("|", "_").replace(":", "_")

    return f"{root}.{index:03d}.{safe_gene}.{safe_region}{ext}"


def main():
    args = parse_args()

    if not os.path.exists(args.input):
        sys.stderr.write(f"[ERROR] Input file not found: {args.input}\n")
        sys.exit(1)

    outdir = os.path.dirname(args.out)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t", dtype=str)

    check_required_columns(df, args.input)

    df = prepare_dataframe(df)

    if df.empty:
        sys.stderr.write("[ERROR] No valid records after filtering.\n")
        sys.exit(1)

    if args.gene:
        df = df[df["Gene"] == args.gene]

    if args.chrom:
        df = df[df["Chrom"] == args.chrom]

    if args.region:
        df = df[df["Region"] == args.region]

    if df.empty:
        sys.stderr.write("[ERROR] No records matched the filtering conditions.\n")
        sys.exit(1)

    group_cols = ["SampleID", "Gene", "Transcript", "Region", "Chrom"]

    groups = list(df.groupby(group_cols, sort=False))

    total = len(groups)

    generated = []

    for i, (key, subdf) in enumerate(groups, start=1):
        _, gene, _, region, _ = key
        output_path = build_output_path(args.out, gene, region, i, total)

        plot_one_group(subdf, args, output_path)
        generated.append(output_path)

    sys.stderr.write("[INFO] Depth ratio plotting finished\n")
    sys.stderr.write(f"[INFO] Input records : {df.shape[0]}\n")
    sys.stderr.write(f"[INFO] Figures       : {len(generated)}\n")

    for p in generated:
        sys.stderr.write(f"[INFO] Output        : {p}\n")


if __name__ == "__main__":
    main()

