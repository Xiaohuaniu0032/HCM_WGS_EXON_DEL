#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot window-level depth ratio for each core gene.

Input:
    *.depth_candidates.all_window_ratio.tsv

Required columns:
    Gene
    Window_ID
    Start
    End
    Depth_Ratio
    Window_Status

Output:
    A multi-page PDF file. Each page shows one gene.

Rules:
    - Window_Status == Normal_window: normal window
    - Window_Status != Normal_window: highlighted in red
    - Depth_Ratio cutoff is fixed at 0.65

Note:
    The upstream workflow already restricts analysis to genes listed in
    HCM_CORE_GENE_LIST, so this script directly plots the supplied input file.
"""

import argparse
import os
import sys

import matplotlib

# Use a non-interactive backend for Linux servers.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages


NORMAL_STATUS = "Normal_window"
RATIO_CUTOFF = 0.65
MAX_LABELS = 30
DPI = 300


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot window-level depth ratio for core genes."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input all_window_ratio.tsv file.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output multi-page PDF file. Each page shows one gene.",
    )

    return parser.parse_args()


def check_required_columns(df):
    required_cols = [
        "Gene",
        "Window_ID",
        "Start",
        "End",
        "Depth_Ratio",
        "Window_Status",
    ]

    missing = [
        col
        for col in required_cols
        if col not in df.columns
    ]

    if missing:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(missing)
            + "\nAvailable columns: "
            + ", ".join(df.columns)
        )


def read_input(input_file):
    if not os.path.isfile(input_file):
        raise FileNotFoundError(
            f"Input file not found: {input_file}"
        )

    if os.path.getsize(input_file) == 0:
        raise ValueError(
            f"Input file is empty: {input_file}"
        )

    df = pd.read_csv(
        input_file,
        sep="\t",
        dtype=str,
    )

    check_required_columns(df)

    df["Start"] = pd.to_numeric(
        df["Start"],
        errors="coerce",
    )

    df["End"] = pd.to_numeric(
        df["End"],
        errors="coerce",
    )

    df["Depth_Ratio"] = pd.to_numeric(
        df["Depth_Ratio"],
        errors="coerce",
    )

    df = df.dropna(
        subset=[
            "Gene",
            "Window_ID",
            "Start",
            "End",
            "Depth_Ratio",
            "Window_Status",
        ]
    ).copy()

    if df.empty:
        raise ValueError(
            "No valid records remain after checking required "
            "fields and numeric values."
        )

    df["Gene"] = (
        df["Gene"]
        .astype(str)
        .str.strip()
    )

    df["Window_ID"] = (
        df["Window_ID"]
        .astype(str)
        .str.strip()
    )

    df["Window_Status"] = (
        df["Window_Status"]
        .astype(str)
        .str.strip()
    )

    df = df[
        (df["Gene"] != "")
        & (df["Window_ID"] != "")
        & (df["Window_Status"] != "")
    ].copy()

    if df.empty:
        raise ValueError(
            "No valid records remain after removing blank fields."
        )

    df = df.sort_values(
        [
            "Gene",
            "Start",
            "End",
            "Window_ID",
        ]
    ).reset_index(drop=True)

    return df


def make_non_normal_text(non_normal_df):
    if non_normal_df.empty:
        return "Non-normal windows: none"

    lines = ["Non-normal windows"]

    for _, row in non_normal_df.head(MAX_LABELS).iterrows():
        lines.append(
            f"{row['Window_ID']}: "
            f"{row['Depth_Ratio']:.3f}"
        )

    if len(non_normal_df) > MAX_LABELS:
        remaining = len(non_normal_df) - MAX_LABELS
        lines.append(f"... +{remaining} more")

    return "\n".join(lines)


def add_annotation_box(ax, text):
    ax.text(
        0.01,
        0.98,
        text,
        transform=ax.transAxes,
        fontsize=7,
        verticalalignment="top",
        horizontalalignment="left",
        bbox=dict(
            boxstyle="round",
            facecolor="white",
            edgecolor="gray",
            alpha=0.88,
        ),
    )


def plot_one_gene(gene_df, gene_name):
    gene_df = gene_df.copy()

    gene_df = gene_df.sort_values(
        [
            "Start",
            "End",
            "Window_ID",
        ]
    ).reset_index(drop=True)

    gene_df["Plot_Index"] = range(
        1,
        len(gene_df) + 1,
    )

    normal_df = gene_df[
        gene_df["Window_Status"] == NORMAL_STATUS
    ]

    non_normal_df = gene_df[
        gene_df["Window_Status"] != NORMAL_STATUS
    ]

    fig_width = max(
        10,
        min(
            24,
            len(gene_df) * 0.12,
        ),
    )

    fig, ax = plt.subplots(
        figsize=(fig_width, 5.8)
    )

    ax.plot(
        gene_df["Plot_Index"],
        gene_df["Depth_Ratio"],
        color="gray",
        linewidth=0.7,
        alpha=0.55,
        zorder=1,
    )

    ax.scatter(
        normal_df["Plot_Index"],
        normal_df["Depth_Ratio"],
        s=18,
        color="lightgray",
        edgecolors="none",
        label="Normal window",
        zorder=2,
    )

    ax.scatter(
        non_normal_df["Plot_Index"],
        non_normal_df["Depth_Ratio"],
        s=36,
        color="red",
        edgecolors="black",
        linewidths=0.3,
        label="Non-normal window",
        zorder=3,
    )

    ax.axhline(
        y=RATIO_CUTOFF,
        color="red",
        linestyle="--",
        linewidth=1,
        alpha=0.75,
        label=(
            f"Depth ratio cutoff = "
            f"{RATIO_CUTOFF}"
        ),
    )

    annotation_text = make_non_normal_text(
        non_normal_df
    )

    add_annotation_box(
        ax,
        annotation_text,
    )

    ax.set_title(
        f"{gene_name}: Window-level Depth Ratio",
        fontsize=14,
        fontweight="bold",
    )

    ax.set_xlabel("Window Order")
    ax.set_ylabel("Depth Ratio")

    ymax = gene_df["Depth_Ratio"].max() + 0.2
    ymax = max(1.2, ymax)

    ax.set_ylim(
        0,
        ymax,
    )

    ax.grid(
        axis="y",
        linestyle="--",
        alpha=0.3,
    )

    ax.legend(
        loc="upper right",
        fontsize=8,
        frameon=True,
    )

    if len(gene_df) <= 60:
        ax.set_xticks(
            gene_df["Plot_Index"]
        )

        ax.set_xticklabels(
            gene_df["Window_ID"],
            rotation=90,
            fontsize=6,
        )
    else:
        ax.tick_params(
            axis="x",
            labelsize=7,
        )

    plt.tight_layout()

    return fig


def write_per_gene_pdf(df, output_pdf):
    output_dir = os.path.dirname(
        os.path.abspath(output_pdf)
    )

    os.makedirs(
        output_dir,
        exist_ok=True,
    )

    gene_count = df["Gene"].nunique()

    print(
        f"[INFO] Start plotting depth ratio "
        f"for {gene_count} core gene(s)."
    )

    with PdfPages(output_pdf) as pdf:
        for gene_name, gene_df in df.groupby(
            "Gene",
            sort=True,
        ):
            fig = plot_one_gene(
                gene_df,
                gene_name,
            )

            pdf.savefig(
                fig,
                dpi=DPI,
            )

            plt.close(fig)

    print(
        f"[INFO] PDF written to: {output_pdf}"
    )


def main():
    args = parse_args()

    try:
        df = read_input(
            args.input
        )

        write_per_gene_pdf(
            df=df,
            output_pdf=args.output,
        )

    except Exception as exc:
        print(
            f"[ERROR] {exc}",
            file=sys.stderr,
        )

        sys.exit(1)


if __name__ == "__main__":
    main()

