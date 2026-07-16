import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

INPUT = ROOT / "tables" / "supplementary_tables" / "Supplementary_Table_S2_Top_50_Motifs.csv"

OUT_PNG = ROOT / "figures" / "main_figures" / "Figure4_FINAL.png"
OUT_TIFF = ROOT / "figures" / "submission_tiff" / "Figure4_Top_Depleted_Motifs.tiff"
OUT_DATA = ROOT / "data" / "submission" / "Figure4_Data.csv"

df = pd.read_csv(INPUT).head(10).copy()

df["gc_content"] = df["kmer"].apply(
    lambda s: (s.count("G") + s.count("C")) / len(s)
)

df["minus_log10_oe"] = -np.log10(df["oe_ratio"])

OUT_DATA.parent.mkdir(parents=True, exist_ok=True)

df.to_csv(OUT_DATA, index=False)

fig = plt.figure(figsize=(12, 5))

ax1 = plt.subplot(1, 2, 1)

colors = plt.cm.Blues(
    (df["gc_content"] - 0.4) / (0.8 - 0.4)
)

bars = ax1.barh(
    df["kmer"],
    df["minus_log10_oe"],
    color=colors
)

ax1.invert_yaxis()

ax1.set_title(
    "A. Ten most severely depleted STR motifs",
    fontweight="bold"
)

ax1.set_xlabel("-log10(Observed/Expected)")

for i, value in enumerate(df["fold_reduction"]):
    ax1.text(
        df["minus_log10_oe"].iloc[i] + 0.05,
        i,
        f"{value:,.0f}",
        va="center",
        fontsize=9
    )

ax2 = plt.subplot(1, 2, 2)

scatter = ax2.scatter(
    df["oe_ratio"],
    df["k"],
    c=df["gc_content"],
    s=220,
    cmap="Blues",
    edgecolor="black"
)

ax2.set_xscale("log")

ax2.set_title(
    "B. Observed/Expected (O/E) ratios of depleted motifs",
    fontweight="bold"
)

ax2.set_xlabel("Observed / Expected (O/E) Ratio")
ax2.set_ylabel("Motif Length (bp)")

for _, row in df.head(3).iterrows():
    ax2.annotate(
        row["kmer"],
        (row["oe_ratio"], row["k"]),
        xytext=(5, 12),
        textcoords="offset points",
        fontsize=8
    )

cbar = plt.colorbar(scatter, ax=ax2)

cbar.set_label("GC Content")

plt.tight_layout()

OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
OUT_TIFF.parent.mkdir(parents=True, exist_ok=True)

plt.savefig(
    OUT_PNG,
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    OUT_TIFF,
    dpi=600,
    bbox_inches="tight"
)

print("Saved:", OUT_PNG)
print("Saved:", OUT_TIFF)
print("Saved:", OUT_DATA)