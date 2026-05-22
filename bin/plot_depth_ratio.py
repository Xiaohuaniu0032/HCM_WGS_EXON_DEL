#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot window-level depth ratio for each gene.

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

Important:
  This script reads the config file first.

  If ANALYZE_CORE_GENES_ONLY=1:
      Plot depth ratio.

  If ANALYZE_CORE_GENES_ONLY=0:
      Skip plotting directly.

Reason:
  In whole-BAM mode, *.all_window_ratio.tsv can be very large,
  and plotting all windows is not practical.
"""

import argparse
import sys

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


NORMAL_STATUS = "Normal_window"
RATIO_CUTOFF = 0.65
MAX_LABELS = 30
DPI = 300


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot gene-level window depth ratio."
    )

    parser.add_argument(
        "--conf",
        required=True,
        help="Config file. ANALYZE_CORE_GENES_ONLY is read from this file."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input all_window_ratio.tsv file."
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output multi-page PDF file. Each page shows one gene."
    )

    return parser.parse_args()


def read_config(conf_file):
    conf = {}

    try:
        with open(conf_file, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip("\n").rstrip("\r")

                if not line.strip():
                    continue

                if line.lstrip().startswith("#"):
                    continue

                # Remove simple inline comments:
                # KEY=value  # comment
                if "#" in line:
                    line = line.split("#", 1)[0].rstrip()

                if "=" not in line:
                    continue

                key, value = line.split("=", 1)

                key = key.strip()
                value = value.strip()

                value = value.strip("'").strip('"')

                if key:
                    conf[key] = value

    except FileNotFoundError:
        raise FileNotFoundError(f"Config file not found: {conf_file}")

    except PermissionError:
        raise PermissionError(f"Permission denied when reading config file: {conf_file}")

    return conf


def normalize_bool(value):
    if value is None:
        return False

    value = str(value).strip().lower()

    if value in {"1", "yes", "true", "on"}:
        return True

    if value in {"0", "no", "false", "off", ""}:
        return False

    return bool(value)


def should_plot(conf_file):
    conf = read_config(conf_file)

    analyze_core_genes_only = normalize_bool(
        conf.get("ANALYZE_CORE_GENES_ONLY", "0")
    )

    if not analyze_core_genes_only:
        print(
            "[INFO] ANALYZE_CORE_GENES_ONLY is not enabled. "
            "Skip plot_depth_ratio.py to avoid plotting whole-BAM depth results."
        )
        return False

    print("[INFO] ANALYZE_CORE_GENES_ONLY=1. Start plotting depth ratio.")
    return True


def check_required_columns(df):
    required_cols = [
        "Gene",
        "Window_ID",
        "Start",
        "End",
        "Depth_Ratio",
        "Window_Status",
    ]

    missing = [col for col in required_cols if col not in df.columns]

    if missing:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(missing)
            + "\nAvailable columns: "
            + ", ".join(df.columns)
        )


def read_input(input_file):
    df = pd.read_csv(input_file, sep="\t", dtype=str)

    check_required_columns(df)

    df["Start"] = pd.to_numeric(df["Start"], errors="coerce")
    df["End"] = pd.to_numeric(df["End"], errors="coerce")
    df["Depth_Ratio"] = pd.to_numeric(df["Depth_Ratio"], errors="coerce")

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

    df = df.sort_values(
        ["Gene", "Start", "End"]
    ).reset_index(drop=True)

    return df


def format_window_status(status):
    """
    Convert internal status to display-friendly text.
    """
    return str(status).replace("_", " ")


def make_non_normal_text(non_normal_df):
    if non_normal_df.empty:
        return "Non-normal windows: none"

    lines = ["Non-normal windows"]

    for _, row in non_normal_df.head(MAX_LABELS).iterrows():
        lines.append(
            f"{row['Window_ID']}: {row['Depth_Ratio']:.3f}"
        )

    if len(non_normal_df) > MAX_LABELS:
        lines.append(f"... +{len(non_normal_df) - MAX_LABELS} more")

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
        ["Start", "End"]
    ).reset_index(drop=True)

    gene_df["Plot_Index"] = range(1, len(gene_df) + 1)

    normal_df = gene_df[gene_df["Window_Status"] == NORMAL_STATUS]
    non_normal_df = gene_df[gene_df["Window_Status"] != NORMAL_STATUS]

    fig_width = max(10, min(24, len(gene_df) * 0.12))

    fig, ax = plt.subplots(figsize=(fig_width, 5.8))

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
        label=f"Depth ratio cutoff = {RATIO_CUTOFF}",
    )

    text = make_non_normal_text(non_normal_df)
    add_annotation_box(ax, text)

    ax.set_title(
        f"{gene_name}: Window-level Depth Ratio",
        fontsize=14,
        fontweight="bold",
    )

    ax.set_xlabel("Window Order")
    ax.set_ylabel("Depth Ratio")

    ymax = gene_df["Depth_Ratio"].max() + 0.2
    ymax = max(1.2, ymax)

    ax.set_ylim(0, ymax)

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
        ax.set_xticks(gene_df["Plot_Index"])
        ax.set_xticklabels(
            gene_df["Window_ID"],
            rotation=90,
            fontsize=6,
        )
    else:
        ax.tick_params(axis="x", labelsize=7)

    plt.tight_layout()

    return fig


def write_per_gene_pdf(df, output_pdf):
    with PdfPages(output_pdf) as pdf:
        for gene_name, gene_df in df.groupby("Gene", sort=True):
            fig = plot_one_gene(gene_df, gene_name)
            pdf.savefig(fig, dpi=DPI)
            plt.close(fig)

    print(f"[INFO] PDF written to: {output_pdf}")


def main():
    args = parse_args()

    if not should_plot(args.conf):
        sys.exit(0)

    df = read_input(args.input)

    write_per_gene_pdf(
        df=df,
        output_pdf=args.output,
    )


if __name__ == "__main__":
    main()

    