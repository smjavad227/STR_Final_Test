import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from config.paths import *

"""
Final combined analysis: T1D (IFNα) + T2D (GSE21232)
"""

import pandas as pd
from scipy.stats import ttest_ind
import numpy as np
import os

DATA = os.path.join(DIABETES_DIR, "data", "methylation")

print("="*60)
print("FINAL COMBINED METHYLATION ANALYSIS")
print("="*60)

# ============================================================
# 1. T1D Analysis (from GSE124809 - IFNα model)
# ============================================================
print("\n1. T1D METHYLATION (IFNα vs CONTROL)")
t1d_df = pd.read_csv(os.path.join(DATA, "methylation_results_t1d_ifna.csv"))

def assign_t1d_group(sample):
    if 'S01' in sample or 'S03' in sample or 'S05' in sample:
        return 'control'
    elif 'S02' in sample or 'S04' in sample or 'S06' in sample:
        return 't1d'
    return 'unknown'

t1d_df['group'] = t1d_df['sample'].apply(assign_t1d_group)
t1d_df = t1d_df[t1d_df['group'] != 'unknown']

t1d_results = []
for motif in t1d_df['motif'].unique():
    control_beta = t1d_df[(t1d_df['motif'] == motif) & (t1d_df['group'] == 'control')]['beta']
    disease_beta = t1d_df[(t1d_df['motif'] == motif) & (t1d_df['group'] == 't1d')]['beta']
    if len(control_beta) > 1 and len(disease_beta) > 1:
        _, p_val = ttest_ind(control_beta, disease_beta, equal_var=False)
        t1d_results.append({
            'condition': 'T1D',
            'motif': motif,
            'mean_control': control_beta.mean(),
            'mean_disease': disease_beta.mean(),
            'fold_change': disease_beta.mean() / control_beta.mean() if control_beta.mean() != 0 else np.nan,
            'p_value': p_val,
            'n_control': len(control_beta),
            'n_disease': len(disease_beta)
        })

t1d_df_results = pd.DataFrame(t1d_results)
print(t1d_df_results.to_string())

# ============================================================
# 2. T2D Analysis (from GSE21232 - 27K)
# ============================================================
print("\n2. T2D METHYLATION (T2D vs CONTROL)")
t2d_df = pd.read_csv(os.path.join(DATA, "methylation_results_t2d.csv"))

control_samples = [f"GSM53088{i}" for i in range(3, 14)]
t2d_samples = [f"GSM53088{i}" for i in range(14, 19)]

t2d_df['group'] = t2d_df['sample'].apply(
    lambda x: 'control' if x in control_samples else 't2d' if x in t2d_samples else 'unknown'
)
t2d_df = t2d_df[t2d_df['group'] != 'unknown']

t2d_results = []
for motif in t2d_df['motif'].unique():
    control_beta = t2d_df[(t2d_df['motif'] == motif) & (t2d_df['group'] == 'control')]['beta']
    disease_beta = t2d_df[(t2d_df['motif'] == motif) & (t2d_df['group'] == 't2d')]['beta']
    if len(control_beta) > 1 and len(disease_beta) > 1:
        _, p_val = ttest_ind(control_beta, disease_beta, equal_var=False)
        t2d_results.append({
            'condition': 'T2D',
            'motif': motif,
            'mean_control': control_beta.mean(),
            'mean_disease': disease_beta.mean(),
            'fold_change': disease_beta.mean() / control_beta.mean() if control_beta.mean() != 0 else np.nan,
            'p_value': p_val,
            'n_control': len(control_beta),
            'n_disease': len(disease_beta)
        })

t2d_df_results = pd.DataFrame(t2d_results)
print(t2d_df_results.to_string())

# ============================================================
# 3. Combine and Save
# ============================================================
combined = pd.concat([t1d_df_results, t2d_df_results], ignore_index=True)
combined.to_csv(os.path.join(DATA, "methylation_combined_stats.csv"), index=False)

print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"T1D motifs: {len(t1d_df_results)}")
print(f"T2D motifs: {len(t2d_df_results)}")
print("\nResults saved to: methylation_combined_stats.csv")
print("="*60)