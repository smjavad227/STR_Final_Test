"""
Figure 05 Generation Script
---------------------------
This script generates Figure 5, visualizing STR motif statistics.
It includes:
1) A bar plot for motif counts by motif length
2) A boxplot with log-scaled observed/expected ratios for motif lengths 4 and 5

Requirements:
    - pandas
    - matplotlib
    - numpy
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def generate_figure_05() -> None:
    """Generate and save Figure 5 in TIFF and PNG formats."""

    # Resolve project paths relative to this script.
    project_root = Path(__file__).resolve().parents[2]
    data_dir = project_root / "figure5_results"
    output_dir = project_root / "figures" / "publication"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Publication-style plotting.
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    # Input files.
    counts_path = data_dir / "figure5_length_stats.csv"
    raw_data_path = data_dir / "figure5_raw_data.csv"

    if not counts_path.exists():
        raise FileNotFoundError(f"Missing file: {counts_path}")
    if not raw_data_path.exists():
        raise FileNotFoundError(f"Missing file: {raw_data_path}")

    # Read data.
    counts_df = pd.read_csv(counts_path)
    raw_df = pd.read_csv(raw_data_path)

    # Expected columns:
    # counts_df: length, count
    # raw_df: length, oe_ratio
    if not {"length", "count"}.issubset(counts_df.columns):
        raise ValueError("figure5_length_stats.csv must contain 'length' and 'count' columns.")
    if not {"length", "oe_ratio"}.issubset(raw_df.columns):
        raise ValueError("figure5_raw_data.csv must contain 'length' and 'oe_ratio' columns.")

    counts_df = counts_df.set_index("length")["count"]

    # Create figure with two panels.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    # Colors: green for motif length 4, red for motif length 5.
    colors = ["#4daf4a", "#e41a1c"]

    # ---------------------------
    # Panel A: Bar plot
    # ---------------------------
    count_4 = int(counts_df.loc[4])
    count_5 = int(counts_df.loc[5])

    bars = ax1.bar(
        ["4", "5"],
        [count_4, count_5],
        color=colors,
        edgecolor="black",
        width=0.6,
        alpha=0.85,
    )

    ax1.set_ylabel("Number of STR motifs")
    ax1.set_xlabel("Motif Length (bp)")
    ax1.set_ylim(0, 1150)
    ax1.set_yticks([0, 200, 400, 600, 800, 1000])
    ax1.grid(axis="y", linestyle="--", alpha=0.5, linewidth=0.8)
    ax1.set_axisbelow(True)

    # Panel label.
    ax1.text(
        0.02,
        1.03,
        "A",
        transform=ax1.transAxes,
        fontsize=13,
        fontweight="bold",
        va="bottom",
        ha="left",
    )

    # Bar labels:
    # 4 bp -> 256
    # 5 bp -> 1280
    ax1.text(
        bars[0].get_x() + bars[0].get_width() / 2,
        count_4 + 20,
        "256",
        ha="center",
        va="bottom",
        fontsize=9,
    )
    ax1.text(
        bars[1].get_x() + bars[1].get_width() / 2,
        count_5 + 20,
        "1280",
        ha="center",
        va="bottom",
        fontsize=9,
    )

    # ---------------------------
    # Panel B: Boxplot + scatter
    # ---------------------------
    data_4 = raw_df.loc[raw_df["length"] == 4, "oe_ratio"].dropna().to_numpy()
    data_5 = raw_df.loc[raw_df["length"] == 5, "oe_ratio"].dropna().to_numpy()

    bp = ax2.boxplot(
        [data_4, data_5],
        positions=[1, 2],
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.5},
        boxprops={"linewidth": 1.2},
        whiskerprops={"linewidth": 1.2},
        capprops={"linewidth": 1.2},
    )

    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    rng = np.random.default_rng(42)
    for i, data in enumerate([data_4, data_5], start=1):
        x = rng.normal(loc=i, scale=0.05, size=len(data))
        ax2.scatter(
            x,
            data,
            color=colors[i - 1],
            alpha=0.2,
            s=10,
            edgecolors="none",
        )

    ax2.set_yscale("log")
    ax2.set_ylim(1e-5, 1)
    ax2.set_ylabel("Observed/Expected Ratio")
    ax2.set_xlabel("Motif Length (bp)")
    ax2.set_xticks([1, 2])
    ax2.set_xticklabels(["4", "5"])
    ax2.grid(axis="y", linestyle="--", alpha=0.5, linewidth=0.8)
    ax2.set_axisbelow(True)

    ax2.text(
        0.02,
        1.03,
        "B",
        transform=ax2.transAxes,
        fontsize=13,
        fontweight="bold",
        va="bottom",
        ha="left",
    )

    # Panel B sample-size annotations.
    ax2.text(1, 2e-5, "n=256", ha="center", fontsize=9)
    ax2.text(2, 2e-5, "n=1024", ha="center", fontsize=9)

    plt.tight_layout(pad=2.0)

    # Save outputs.
    tiff_path = output_dir / "Figure05.tiff"
    png_path = output_dir / "Figure05.png"

    fig.savefig(
        tiff_path,
        dpi=300,
        format="tiff",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    fig.savefig(png_path, dpi=300)

    plt.close(fig)

    print("Figure 5 saved successfully:")
    print(f"  TIFF: {tiff_path}")
    print(f"  PNG : {png_path}")


if __name__ == "__main__":
    generate_figure_05()
