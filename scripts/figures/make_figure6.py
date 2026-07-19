import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from config.paths import *

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_FILE = PROJECT_ROOT / "data" / "figure6" / "Figure6_Raw_Data.csv"
OUT_DIR = PROJECT_ROOT / "figures"

OUT_DIR.mkdir(exist_ok=True)

df = pd.read_csv(DATA_FILE)

metric_labels = {
    "length": "Length\n(bp)",
    "period": "Period\n(bp)",
    "copies": "Copy\nNumber",
    "score": "TRF\nScore",
    "entropy": "Entropy\n(bits)"
}

labels = [metric_labels[m] for m in df["metric"]]

# --------------------------------------------------
# Figure 6A
# --------------------------------------------------

x = np.arange(len(df))
width = 0.36

fig, ax = plt.subplots(figsize=(10,6))

bars1 = ax.bar(
    x - width/2,
    df["x_mean"],
    width,
    label="Chromosome X"
)

bars2 = ax.bar(
    x + width/2,
    df["auto_mean"],
    width,
    label="Autosomes (chr1,8,19,21)"
)

ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel("Mean Value")
ax.set_title(
    "A. Comparison of STR Characteristics: Chromosome X vs Autosomes"
)

ax.legend()
ax.grid(axis="y", linestyle="--", alpha=0.3)

for bars in [bars1, bars2]:
    for bar in bars:
        h = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width()/2,
            h,
            f"{h:.1f}",
            ha="center",
            va="bottom",
            fontsize=9
        )

plt.tight_layout()

plt.savefig(
    OUT_DIR / "Figure6A.png",
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    OUT_DIR / "Figure6A.tiff",
    dpi=600,
    bbox_inches="tight"
)

plt.close()

# --------------------------------------------------
# Figure 6B
# --------------------------------------------------

fig, ax = plt.subplots(figsize=(8, 7))

xvals = df["x_mean"].values
yvals = df["auto_mean"].values

labels_plot = ["Length", "Period", "Copies", "Score", "Entropy"]
markers = ["o", "s", "^", "D", "v"]
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]

for i in range(len(df)):
    ax.scatter(
        xvals[i],
        yvals[i],
        s=400,
        marker=markers[i],
        color=colors[i],
        alpha=0.8,
        label=labels_plot[i]
    )

slope, intercept, r, p, se = stats.linregress(xvals, yvals)

xmin = 0
xmax = max(xvals.max(), yvals.max()) * 1.2

line_x = np.linspace(xmin, xmax, 200)

# y=x
ax.plot(
    line_x,
    line_x,
    "--",
    linewidth=2,
    alpha=0.6,
    label="y = x"
)

# regression
ax.plot(
    line_x,
    slope * line_x + intercept,
    linewidth=3,
    alpha=0.9,
    label=f"Regression (R² = {r*r:.3f})"
)

# stats box
ax.text(
    0.05,
    0.95,
    f"Pearson r = {r:.3f}\np = {p:.4f}",
    transform=ax.transAxes,
    fontsize=13,
    verticalalignment="top",
    bbox=dict(
        boxstyle="round",
        facecolor="white",
        edgecolor="black"
    )
)

ax.set_xlabel("Chromosome X (mean value)", fontsize=14)
ax.set_ylabel("Autosomes (mean value)", fontsize=14)

ax.set_title(
    "B. Correlation between Chromosome X and Autosomes",
    fontsize=18,
    fontweight="bold",
    pad=15
)

ax.grid(True, linestyle="--", alpha=0.3)

ax.legend(
    loc="lower right",
    fontsize=11
)

plt.tight_layout()

plt.savefig(
    OUT_DIR / "Figure6B.png",
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    OUT_DIR / "Figure6B.tiff",
    dpi=600,
    bbox_inches="tight"
)

plt.close()