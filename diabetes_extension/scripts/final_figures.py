"""
final_figures.py

Generates both final figures:
- Figure 3: Heatmap of methylation fold changes (T1D and T2D)
- Figure 4: Bar plot of expression and methylation fold changes with significance stars
Uses the final confirmed data files.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

BASE = r"D:\STR_Diabetes_Extension"
DATA_DIR = os.path.join(BASE, "data", "methylation")
FIGURES_DIR = os.path.join(BASE, "results", "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

# ============================================================================
# 1. Load methylation final stats
# ============================================================================
methyl_stats = pd.read_csv(os.path.join(DATA_DIR, "methylation_final_stats.csv"))
methyl_stats = methyl_stats[['condition', 'motif', 'fold_change', 'p_value']].drop_duplicates()
methyl_stats['condition'] = methyl_stats['condition'].replace({'T1D': 'T1D', 'T2D': 'T2D'})

# ============================================================================
# 2. FIGURE 3: Heatmap
# ============================================================================
pivot = methyl_stats.pivot(index='motif', columns='condition', values='fold_change')
if 'T2D' in pivot.columns:
    pivot = pivot.sort_values('T2D', ascending=True)

plt.figure(figsize=(6, 8))
sns.heatmap(pivot, annot=True, fmt='.2f', cmap='YlOrRd',
            cbar_kws={'label': 'Fold Change (Disease / Control)'},
            linewidths=0.5, linecolor='gray')
plt.title('Methylation Fold Changes at Conserved Non-CG Motifs', fontsize=14)
plt.ylabel('Motif')
plt.xlabel('Disease')
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "Figure3_Methylation_FoldChange_Heatmap.png"), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIGURES_DIR, "Figure3_Methylation_FoldChange_Heatmap.tiff"), dpi=300, bbox_inches='tight', format='tiff')
plt.close()

# ============================================================================
# 3. FIGURE 4: Bar plot with stars
# ============================================================================
expression_data = {
    'Condition': ['Normal', 'T1D', 'T2D'],
    'Fold_Change': [1.53, 1.11, 0.84],
    'p_value': [0.039, 1.00, 1.00],
    'Type': ['Expression', 'Expression', 'Expression']
}

avg_methyl = methyl_stats.groupby('condition').agg({
    'fold_change': 'mean',
    'p_value': lambda x: x.min()
}).reset_index()
avg_methyl.columns = ['Condition', 'Fold_Change', 'p_value']
avg_methyl['Type'] = 'Methylation'
avg_methyl['Condition'] = avg_methyl['Condition'].replace({'T1D': 'T1D (Methyl)', 'T2D': 'T2D (Methyl)'})

bar_df = pd.concat([pd.DataFrame(expression_data), avg_methyl], ignore_index=True)

plt.figure(figsize=(10, 6))
sns.barplot(data=bar_df, x='Condition', y='Fold_Change', hue='Type', palette='Set2')
plt.axhline(y=1, color='gray', linestyle='--', label='No Change')

for i, row in bar_df.iterrows():
    if row['p_value'] < 0.05:
        if row['p_value'] < 0.001:
            star = '***'
        elif row['p_value'] < 0.01:
            star = '**'
        else:
            star = '*'
        plt.text(i, row['Fold_Change'] + 0.05, star, ha='center', fontsize=14, fontweight='bold')

plt.title('Fold Changes: Expression and Methylation', fontsize=14)
plt.ylabel('Fold Change (Disease / Control)')
plt.xlabel('')
plt.legend(title='Data Type')
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "Figure4_FoldChange_Barplot.png"), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIGURES_DIR, "Figure4_FoldChange_Barplot.tiff"), dpi=300, bbox_inches='tight', format='tiff')
plt.close()

# ============================================================================
# 4. Done
# ============================================================================
print("Both figures regenerated successfully.")
print("Figure 3: ", os.path.join(FIGURES_DIR, "Figure3_Methylation_FoldChange_Heatmap.png"))
print("Figure 4: ", os.path.join(FIGURES_DIR, "Figure4_FoldChange_Barplot.png"))