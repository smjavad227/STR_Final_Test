import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.paths import *

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Figure 3
--------
Heatmap of average conservation across primate genomes for the
top 50 STR motifs ranked by average conservation.

Input
-----
Figure3_Paper_Submission/Figure3_FullData.csv

Outputs
-------
Figure3.tiff
Figure3.pdf
Figure3_preview.png
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

FIG_WIDTH = 3.35
FIG_HEIGHT = 4.00
DPI = 600

TOP_N = 50

HIGHLIGHT_MOTIFS = {
    "ACAC",
    "AAGA",
    "TTTAA",
    "TCTT",
    "GAGG",
    "TATC",
    "ATCT",
    "ATAG",
    "TATG",
    "AGAT",
    "AAATT",
}

SPECIES_COLUMNS = [
    "Human",
    "Chimpanzee",
    "Gorilla",
    "Orangutan",
    "Macaque",
]


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent

INPUT_FILE = (
    PROJECT_ROOT
    / "Figure3_Paper_Submission"
    / "Figure3_FullData.csv"
)

OUTPUT_TIFF = PROJECT_ROOT / "Figure3_Paper_Submission" / "Figure3.tiff"
OUTPUT_PDF = PROJECT_ROOT / "Figure3_Paper_Submission" / "Figure3.pdf"
OUTPUT_PNG = PROJECT_ROOT / "Figure3_Paper_Submission" / "Figure3_preview.png"


# ---------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------

def load_data(csv_file: Path) -> pd.DataFrame:
    """
    Load the conservation table.
    """

    df = pd.read_csv(csv_file)

    required_columns = [
        "kmer",
        "avg_conservation",
        *SPECIES_COLUMNS,
    ]

    missing = [c for c in required_columns if c not in df.columns]

    if missing:
        raise ValueError(
            f"Missing required columns: {', '.join(missing)}"
        )

    return df


# ---------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------

def prepare_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Select the top motifs ranked by average conservation.
    """

    plot_df = (
        df.sort_values(
            by="avg_conservation",
            ascending=False,
        )
        .head(TOP_N)
        .reset_index(drop=True)
    )

    return plot_df


# ---------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------

def draw_heatmap(ax, plot_df: pd.DataFrame):
    """
    Draw heatmap.
    """

    matrix = plot_df[SPECIES_COLUMNS].to_numpy()

    image = ax.imshow(
        matrix,
        aspect="auto",
        cmap="viridis",
        interpolation="nearest",
        origin="upper",
        vmin=matrix.min(),
        vmax=matrix.max(),
    )

    ax.set_xticks(range(len(SPECIES_COLUMNS)))
    ax.set_xticklabels(
        SPECIES_COLUMNS,
        rotation=45,
        ha="right",
        fontsize=7,
    )

    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(
        plot_df["kmer"],
        fontsize=6,
    )

    ax.tick_params(length=0)

    for spine in ax.spines.values():
        spine.set_visible(False)

    return image


def highlight_selected_rows(ax, plot_df: pd.DataFrame):
    """
    Highlight selected motifs.
    """

    yticklabels = ax.get_yticklabels()

    for row_index, motif in enumerate(plot_df["kmer"]):

        if motif not in HIGHLIGHT_MOTIFS:
            continue

        rectangle = patches.Rectangle(
            (-0.5, row_index - 0.5),
            len(SPECIES_COLUMNS),
            1,
            fill=False,
            edgecolor="black",
            linewidth=1.2,
        )

        ax.add_patch(rectangle)

        yticklabels[row_index].set_fontweight("bold")


def add_colorbar(fig, ax, image):
    """
    Add colorbar.
    """

    cbar = fig.colorbar(
        image,
        ax=ax,
        fraction=0.046,
        pad=0.04,
    )

    cbar.set_label(
        "Conservation",
        fontsize=7,
    )

    cbar.ax.tick_params(
        labelsize=6,
    )


# ---------------------------------------------------------------------
# Figure creation
# ---------------------------------------------------------------------

def create_figure(plot_df: pd.DataFrame):

    fig, ax = plt.subplots(
        figsize=(FIG_WIDTH, FIG_HEIGHT),
    )

    image = draw_heatmap(
        ax=ax,
        plot_df=plot_df,
    )

    highlight_selected_rows(
        ax=ax,
        plot_df=plot_df,
    )

    add_colorbar(
        fig=fig,
        ax=ax,
        image=image,
    )

    plt.tight_layout()

    return fig


# ---------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------

def save_figure(fig):

    fig.savefig(
        OUTPUT_TIFF,
        dpi=DPI,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )

    fig.savefig(
        OUTPUT_PDF,
        dpi=DPI,
        bbox_inches="tight",
    )

    fig.savefig(
        OUTPUT_PNG,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():

    dataframe = load_data(INPUT_FILE)

    plot_dataframe = prepare_data(dataframe)

    figure = create_figure(plot_dataframe)

    save_figure(figure)

    print(f"Input : {INPUT_FILE}")
    print(f"Saved : {OUTPUT_TIFF}")
    print(f"Saved : {OUTPUT_PDF}")
    print(f"Saved : {OUTPUT_PNG}")


if __name__ == "__main__":
    main()