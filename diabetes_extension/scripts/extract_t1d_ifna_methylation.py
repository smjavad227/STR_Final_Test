import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from config.paths import *

"""
Extract methylation data for T1D from GSE124809 (IFNa-treated human islets).
"""

import pandas as pd
import os
from intervaltree import IntervalTree

DATA = os.path.join(DIABETES_DIR, "data", "methylation")

MOTIF_BED = os.path.join(DIABETES_DIR, "data", "target_motifs_filtered.bed")
PROBE_BED = os.path.join(DATA, "450k_probes_hg19.bed")
MATRIX_FILE = os.path.join(DATA, "GSE124809_AVG_Beta_DectectionP_data.txt.gz")
OUTPUT_FILE = os.path.join(DATA, "methylation_results_t1d_ifna.csv")

print("="*60)
print("EXTRACTING T1D METHYLATION DATA (GSE124809)")
print("="*60)

# 1. Load motifs
print("Loading motifs...")
motifs = pd.read_csv(MOTIF_BED, sep='\t', header=None,
                     names=['chr', 'start', 'end', 'motif', 'gene'])
motifs = motifs[motifs['chr'] == 'X']
print(f"Loaded {len(motifs)} motifs on chromosome X.")

# 2. Load 450K probes
print("Loading 450K probes...")
probes = pd.read_csv(PROBE_BED, sep='\t', header=None,
                     names=['chr', 'start', 'end', 'probe_id'])
probes = probes[probes['chr'] == 'X']
print(f"Loaded {len(probes)} probes on chromosome X.")

# 3. Compute overlaps
print("Computing overlaps...")
tree = IntervalTree()
for _, row in probes.iterrows():
    tree.addi(row['start'], row['end'], row['probe_id'])

overlaps = []
for _, motif in motifs.iterrows():
    for interval in tree.overlap(motif['start'], motif['end']):
        overlaps.append({
            'probe_id': interval.data,
            'motif': motif['motif'],
            'gene': motif['gene']
        })

overlap_df = pd.DataFrame(overlaps).drop_duplicates()
print(f"Found {len(overlap_df)} probe-motif overlaps.")

if overlap_df.empty:
    print("No overlaps. Exiting.")
    exit()

# 4. Read GSE124809 matrix
print("Reading GSE124809 matrix...")
df = pd.read_csv(MATRIX_FILE, sep='\t', header=0, compression='gzip')

# Find the probe ID column
probe_col = None
for col in df.columns:
    if 'probe' in col.lower() or 'id' in col.lower():
        probe_col = col
        break

if probe_col is None:
    print("ERROR: Could not find probe ID column.")
    print("Columns:", df.columns.tolist())
    exit()

df.set_index(probe_col, inplace=True)

# 5. Find overlapping probes in matrix
probe_ids = set(overlap_df['probe_id'])
available = [p for p in probe_ids if p in df.index]
print(f"Available probes in matrix: {len(available)}")

if not available:
    print("No overlapping probes found in GSE124809 matrix.")
    exit()

# 6. Extract and melt
sub = df.loc[available]
melted = sub.reset_index().melt(id_vars=probe_col, var_name='sample', value_name='beta')
melted.rename(columns={probe_col: 'probe_id'}, inplace=True)

# 7. Merge with overlap info
result = pd.merge(melted, overlap_df, on='probe_id', how='inner')

# 8. Aggregate by motif and sample
grouped = result.groupby(['sample', 'motif', 'gene']).agg({'beta': 'mean'}).reset_index()
grouped['dataset'] = 'T1D_IFNa'

# 9. Save
grouped.to_csv(OUTPUT_FILE, index=False)
print(f"\nSUCCESS: {len(grouped)} records saved to {OUTPUT_FILE}")
print("Unique motifs:", grouped['motif'].unique())