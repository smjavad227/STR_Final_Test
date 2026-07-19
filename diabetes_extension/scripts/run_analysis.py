import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from config.paths import *

"""
Full differential expression analysis for STR Diabetes Extension project.
Compares target genes (near conserved non-CG motifs) vs control genes
in pancreatic beta cells under Normal, T1D, and T2D conditions.
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')


# =============================================================================
# 1. Path configuration
# =============================================================================
DATA_FILE = os.path.join(DIABETES_DIR_DIR, "data", "expression_data.csv")
OUT_TABLES = os.path.join(DIABETES_DIR_DIR, "results", "tables")
OUT_FIGURES = os.path.join(DIABETES_DIR_DIR, "results", "figures")

os.makedirs(OUT_TABLES, exist_ok=True)
os.makedirs(OUT_FIGURES, exist_ok=True)

print("Directories are ready.")


# =============================================================================
# 2. Data loading and cleaning
# =============================================================================
df = pd.read_csv(DATA_FILE)

# Remove rows where Expression is missing (originally "No data")
df_clean = df.dropna(subset=['Expression']).copy()
print(f"Loaded {len(df)} records, kept {len(df_clean)} after removing missing values.")


# =============================================================================
# 3. Helper function for Cohen's d (effect size)
# =============================================================================
def cohen_d(x, y):
    """
    Calculate Cohen's d for two independent samples.
    Uses pooled standard deviation (equal variance assumption).
    """
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return np.nan
    var_x = np.var(x, ddof=1)
    var_y = np.var(y, ddof=1)
    pooled_std = np.sqrt(((nx - 1) * var_x + (ny - 1) * var_y) / (nx + ny - 2))
    if pooled_std == 0:
        return np.nan
    return (np.mean(x) - np.mean(y)) / pooled_std


# =============================================================================
# 4. Statistical analysis per condition
# =============================================================================
conditions = ['Normal', 'T1D', 'T2D']
results = []

for cond in conditions:
    sub = df_clean[df_clean['Condition'] == cond]
    target_expr = sub[sub['Group'] == 'Target']['Expression']
    control_expr = sub[sub['Group'] == 'Control']['Expression']
    
    n_target = len(target_expr)
    n_control = len(control_expr)
    mean_target = target_expr.mean()
    mean_control = control_expr.mean()
    fold_change = mean_target / mean_control if mean_control != 0 else np.nan
    
    # Statistical tests (only if both groups have at least 2 observations)
    if n_target >= 2 and n_control >= 2:
        t_stat, p_t = ttest_ind(target_expr, control_expr, equal_var=False)
        u_stat, p_mw = mannwhitneyu(target_expr, control_expr, alternative='two-sided')
        d = cohen_d(target_expr, control_expr)
    else:
        p_t, p_mw, d = np.nan, np.nan, np.nan
    
    results.append({
        'Condition': cond,
        'Target_N': n_target,
        'Control_N': n_control,
        'Target_Mean': mean_target,
        'Control_Mean': mean_control,
        'Fold_Change': fold_change,
        'Cohen_d': d,
        'p_value_t_test': p_t,
        'p_value_MW': p_mw
    })

stats_df = pd.DataFrame(results)

# Bonferroni correction for 3 comparisons
stats_df['p_t_adj'] = multipletests(stats_df['p_value_t_test'], method='bonferroni')[1]
stats_df['p_mw_adj'] = multipletests(stats_df['p_value_MW'], method='bonferroni')[1]

# Save table
stats_path = os.path.join(OUT_TABLES, "expression_stats.csv")
stats_df.to_csv(stats_path, index=False)
print(f"Statistical table saved to: {stats_path}")

# Print summary
print("\n" + "="*70)
print("Statistical summary:")
print(stats_df.to_string(float_format="%.4f"))
print("="*70 + "\n")


# =============================================================================
# 5. Figure 1: Boxplot comparing Target vs Control across conditions
# =============================================================================
plt.figure(figsize=(10, 6))
sns.boxplot(data=df_clean, x='Condition', y='Expression', hue='Group', palette='Set2')
plt.title('Target vs Control Gene Expression in Pancreatic Beta Cells', fontsize=14)
plt.ylabel('Scaled Expression (0-1)')
plt.xlabel('Condition')
plt.legend(title='Group')
plt.tight_layout()

fig1_png = os.path.join(OUT_FIGURES, "Figure1_Boxplot.png")
fig1_tiff = os.path.join(OUT_FIGURES, "Figure1_Boxplot.tiff")
plt.savefig(fig1_png, dpi=300, bbox_inches='tight')
plt.savefig(fig1_tiff, dpi=300, bbox_inches='tight', format='tiff')
plt.close()
print("Figure 1 (Boxplot) saved as PNG and TIFF.")


# =============================================================================
# 6. Figure 2: Heatmap of mean expression per gene across conditions
# =============================================================================
# Pivot table: genes as rows, conditions as columns, values = mean expression
pivot = df_clean.pivot_table(index='Gene', columns='Condition', values='Expression', aggfunc='mean')
# Sort genes by Normal expression (high to low) for better visualisation
pivot = pivot.reindex(pivot['Normal'].sort_values(ascending=False).index)

plt.figure(figsize=(8, 10))
sns.heatmap(pivot, annot=True, cmap='viridis', fmt='.2f',
            cbar_kws={'label': 'Mean Scaled Expression'})
plt.title('Gene Expression Heatmap in Beta Cells\n(Normal, T1D, T2D)', fontsize=14)
plt.tight_layout()

fig2_png = os.path.join(OUT_FIGURES, "Figure2_Heatmap.png")
fig2_tiff = os.path.join(OUT_FIGURES, "Figure2_Heatmap.tiff")
plt.savefig(fig2_png, dpi=300, bbox_inches='tight')
plt.savefig(fig2_tiff, dpi=300, bbox_inches='tight', format='tiff')
plt.close()
print("Figure 2 (Heatmap) saved as PNG and TIFF.")


# =============================================================================
# 7. Completion message
# =============================================================================
print("\nAnalysis completed successfully.")
print(f"All outputs are in: {DIABETES_DIR_DIR}")
print(f"  - Table : {OUT_TABLES}")
print(f"  - Figures: {OUT_FIGURES}")