import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.paths import *

import sys
from pathlib import Path

DIABETES_DIR = Path(__file__).resolve().parents[2]

if str(DIABETES_DIR) not in sys.path:
    sys.path.insert(0, str(DIABETES_DIR))

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

import sys; sys.path.append('..'); from config.output_paths import FIG1_DIR

DIABETES_DIR = Path(__file__).resolve().parents[2]

INPUT_FILE = DIABETES_DIR / "results" / "FINAL_evolutionary_forbidden_integration.tsv"

FIG1_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT_FILE, sep="\t")

cpg = df[df["has_CG"] == True]["oe_ratio"]
non_cpg = df[df["has_CG"] == False]["oe_ratio"]

u_stat, p_value = mannwhitneyu(cpg, non_cpg, alternative="two-sided")

fig, ax = plt.subplots(figsize=(6, 6))

data = [cpg, non_cpg]

bp = ax.boxplot(data, patch_artist=True, showfliers=False)

bp["boxes"][0].set_facecolor("#d95f02")
bp["boxes"][1].set_facecolor("#1f77b4")

for median in bp["medians"]:
    median.set_linewidth(2)

rng = np.random.default_rng(42)

for i, values in enumerate(data, start=1):
    x = rng.normal(i, 0.06, len(values))
    ax.scatter(x, values, alpha=0.45, s=18)

means = [cpg.mean(), non_cpg.mean()]

ax.scatter([1, 2], means, marker="D", s=40, edgecolor="black", zorder=5)

ax.set_yscale("log")

ax.set_ylabel("Observed-to-Expected (O/E) Ratio")

ax.set_xticks([1, 2])

ax.set_xticklabels([
    f"CpG-containing\n\n(n={len(cpg)})",
    f"Non-CpG\n\n(n={len(non_cpg)})"
])

ax.set_title(
    "CpG-Containing STR Motifs Exhibit Lower O/E Ratios",
    pad=12,
    weight="bold"
)

ax.text(1.18, 0.07, "Mann–Whitney U test", fontsize=9)
ax.text(1.18, 0.045, f"p = {p_value:.2e}", fontsize=9)

ax.grid(True, alpha=0.3)

plt.tight_layout()

outfile = FIG1_DIR / "Figure1_CpG_Depletion.tiff"

plt.savefig(outfile, dpi=300, format="tiff", pil_kwargs={"compression": "tiff_lzw"})
plt.close()

print(f"Saved: {outfile}")
print(f"p-value = {p_value:.3e}")