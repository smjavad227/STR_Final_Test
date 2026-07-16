"""
Extract 27K probe coordinates from GPL8490 full table (saved from GEO).
Column names: 'Chr' and 'MapInfo' (not 'CHR' and 'MAPINFO').
"""

import pandas as pd
import os
from intervaltree import IntervalTree

# ============================================================================
# Paths
# ============================================================================
BASE_DIR = r"D:\STR_Diabetes_Extension"
DATA_DIR = os.path.join(BASE_DIR, "data", "methylation")

TABLE_FILE = os.path.join(DATA_DIR, "GPL8490_full_table.txt")
PROBE_BED_FILE = os.path.join(DATA_DIR, "27k_probes_hg19.bed")
MOTIF_BED_FILE = os.path.join(BASE_DIR, "data", "target_motifs_filtered.bed")
MATRIX_FILE = os.path.join(DATA_DIR, "GSE21232_series_matrix.txt.gz")
OUTPUT_CSV = os.path.join(DATA_DIR, "methylation_results_t2d.csv")

# ============================================================================
# 1. Read the full table (tab-separated, first row is header)
# ============================================================================
if not os.path.exists(TABLE_FILE):
    print(f"ERROR: Table file not found: {TABLE_FILE}")
    print("Please save the full table from GPL8490 page as GPL8490_full_table.txt")
    exit()

print("Reading GPL8490 full table...")
df = pd.read_csv(TABLE_FILE, sep='\t', comment='#', header=0)

# Use correct column names: 'Chr' and 'MapInfo'
if 'ID' not in df.columns or 'Chr' not in df.columns or 'MapInfo' not in df.columns:
    print("ERROR: Required columns (ID, Chr, MapInfo) not found.")
    print("Available columns:", df.columns.tolist())
    exit()

# ============================================================================
# 2. Extract chromosome X probes
# ============================================================================
print("Extracting chromosome X probes...")
bed = df[['ID', 'Chr', 'MapInfo']].copy()
bed['end'] = bed['MapInfo'] + 1
bed = bed[['Chr', 'MapInfo', 'end', 'ID']]
bed.columns = ['chr', 'start', 'end', 'probe_id']

# Clean chromosome names and filter for X
bed['chr'] = bed['chr'].astype(str).str.replace('chr', '').str.replace('X', 'X')
bed = bed[bed['chr'] == 'X']
bed = bed.dropna()

bed.to_csv(PROBE_BED_FILE, sep='\t', index=False, header=False)
print(f"SUCCESS: Extracted {len(bed)} probes on chromosome X.")
print(f"Saved to: {PROBE_BED_FILE}")

# ============================================================================
# 3. Continue with overlap and methylation analysis
# ============================================================================
print("Loading motifs...")
motifs = pd.read_csv(MOTIF_BED_FILE, sep='\t', header=None,
                     names=['chr', 'start', 'end', 'motif', 'gene'])
motifs = motifs[motifs['chr'] == 'X']
print(f"Loaded {len(motifs)} motifs on chromosome X.")

print("Loading 27K probes...")
probes = pd.read_csv(PROBE_BED_FILE, sep='\t', header=None,
                     names=['chr', 'start', 'end', 'probe_id'])
probes = probes[probes['chr'] == 'X']
print(f"Loaded {len(probes)} probes on chromosome X.")

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
            'gene': motif['gene'],
            'motif_start': motif['start'],
            'motif_end': motif['end']
        })

overlap_df = pd.DataFrame(overlaps).drop_duplicates()
print(f"Found {len(overlap_df)} probe-motif overlaps.")

if overlap_df.empty:
    print("No overlaps. Exiting. (27K platform lacks these probes.)")
    exit()

# Read GSE21232 matrix
if not os.path.exists(MATRIX_FILE):
    print(f"ERROR: Matrix file not found: {MATRIX_FILE}")
    print("Please download GSE21232_series_matrix.txt.gz from GEO.")
    exit()

print("Reading GSE21232 matrix...")
df_matrix = pd.read_csv(MATRIX_FILE, sep='\t', comment='!', header=0)
df_matrix.set_index('ID_REF', inplace=True)

probe_ids = set(overlap_df['probe_id'])
available = [p for p in probe_ids if p in df_matrix.index]
print(f"Available probes in matrix: {len(available)}")

if not available:
    print("No overlapping probes found in GSE21232 matrix.")
    exit()

sub = df_matrix.loc[available]
melted = sub.reset_index().melt(id_vars='ID_REF', var_name='sample', value_name='beta')
melted.rename(columns={'ID_REF': 'probe_id'}, inplace=True)
result = pd.merge(melted, overlap_df, on='probe_id', how='inner')

result.to_csv(OUTPUT_CSV, index=False)
print(f"SUCCESS: {len(result)} records saved to {OUTPUT_CSV}")
print("Unique motifs represented:", result['motif'].unique())