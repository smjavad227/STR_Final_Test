import sys
from pathlib import Path

# =========================
# PROJECT ROOT FIX
# =========================
BASE_DIR = Path(__file__).resolve().parents[2]

if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

# =========================
# IMPORTS
# =========================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, linregress

import sys; sys.path.append('..'); from config.output_paths import FIG2_DIR

# =========================
# INPUT
# =========================
INPUT_FILE = BASE_DIR / "results" / "FINAL_evolutionary_forbidden_integration.tsv"

FIG2_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT_FILE, sep="\t")
df = df[df["conservation_score"].notna()].copy()

r, p = pearsonr(df["oe_ratio"], df["conservation_score"])
slope, intercept, *_ = linregress(df["oe_ratio"], df["conservation_score"])

fig, ax = plt.subplots(figsize=(11, 7))
gs = fig.add_gridspec(1, 2, width_ratios=[4, 1.3], wspace=0.05)

ax = fig.add_subplot(gs[0])
ax_hist = fig.add_subplot(gs[1], sharey=ax)

color_map = {
    "Moderate": "#4C78A8",
    "Strong": "#F28E2B",
    "Very strong": "#E15759",
    "Extreme": "#B279A2"
}

for level in color_map:
    subset = df[df["significance_level"] == level]
    if len(subset) == 0:
        continue

    ax.scatter(
        subset["oe_ratio"],
        subset["conservation_score"],
        s=45,
        alpha=0.8,
        edgecolor="black",
        linewidth=0.4,
        color=color_map[level],
        label=level
    )

xline = np.linspace(df["oe_ratio"].min(), df["oe_ratio"].max(), 200)
yline = slope * xline + intercept
ax.plot(xline, yline, color="black", linewidth=2)

ax.set_title(
    "Association Between STR Motif O/E Ratios and Conservation",
    fontsize=14,
    weight="bold"
)

ax.set_xlabel("O/E Ratio", fontsize=11)
ax.set_ylabel("Conservation Score", fontsize=11)

ax.text(
    0.03, 0.95,
    f"r={r:.3f}\nR²={r**2:.3f}\np={p:.3g}\nn={len(df)}",
    transform=ax.transAxes,
    fontsize=10,
    weight="bold"
)

ax.grid(alpha=0.25)
ax.legend()

ax_hist.hist(
    df["conservation_score"],
    bins=16,
    orientation="horizontal",
    color="#6BAED6",
    edgecolor="black",
    alpha=0.9
)

ax_hist.axhline(df["conservation_score"].mean(), color="green", linestyle="--")
ax_hist.axhline(df["conservation_score"].median(), color="red", linestyle="--")

ax_hist.set_xlabel("Frequency")
plt.setp(ax_hist.get_yticklabels(), visible=False)

plt.tight_layout()

outfile = FIG2_DIR / "Figure2_Conservation_Correlation.tiff"
plt.savefig(outfile, dpi=300, format="tiff")
plt.close()

print("Saved:", outfile)