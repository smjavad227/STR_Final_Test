#!/usr/bin/env bash
# =============================================================================
# 02_annotate_motifs.sh
#
# Annotates the genomic coordinates of eleven conserved non-CG STR motifs
# with chromosomal cytoband locations based on the human genome assembly
# GRCh38/hg38.
#
# Input:
#   data/intermediate/all_11_motifs_coordinates.bed
#     Four-column file with header: #chr, start, end, motif.
#
# Reference:
#   data/reference/cytoBand.txt
#     UCSC cytoBand track for GRCh38/hg38. Downloaded automatically if absent.
#
# Output:
#   data/output/motifs_11_paper_final_all.tsv
#     Five-column table with header:
#     Motif, Chromosome, Start, End, Cytoband
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

INPUT="data/intermediate/all_11_motifs_coordinates.bed"
CYTO_TXT="data/reference/cytoBand.txt"
CYTO_GZ="data/intermediate/cytoBand.txt.gz"
NOHEADER="data/intermediate/motifs_noheader.bed"
STD_CHR="data/intermediate/motifs_standardized_chr.bed"
INTERSECTED="data/intermediate/motifs_cytoband_intersected.bed"
FINAL="data/output/motifs_11_paper_final_all.tsv"

mkdir -p data/intermediate data/output data/reference

# Validate input
if [[ ! -s "${INPUT}" ]]; then
    echo "ERROR: input file not found or empty: ${INPUT}"
    exit 1
fi

# -----------------------------------------------------------------------------
# Step 1: Remove header from input
# -----------------------------------------------------------------------------
echo "Step 1: Removing header from input"
tail -n +2 "${INPUT}" > "${NOHEADER}"
echo "  Records: $(wc -l < "${NOHEADER}")"

# -----------------------------------------------------------------------------
# Step 2: Download and prepare cytoBand reference if absent
# -----------------------------------------------------------------------------
if [[ ! -s "${CYTO_TXT}" ]]; then
    echo "Step 2: Downloading UCSC cytoBand reference for hg38"
    wget -q -O "${CYTO_GZ}" \
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"
    python3 - <<PY
import gzip, shutil
with gzip.open("${CYTO_GZ}", "rb") as fin, open("${CYTO_TXT}", "wb") as fout:
    shutil.copyfileobj(fin, fout)
PY
    echo "  Extracted: ${CYTO_TXT}"
else
    echo "Step 2: Using existing reference: ${CYTO_TXT}"
fi

# -----------------------------------------------------------------------------
# Step 3: Standardize chromosome nomenclature (X -> chrX)
# The UCSC reference uses the chr prefix, whereas the input uses X.
# Standardizing avoids missed intersections for chromosome X records.
# -----------------------------------------------------------------------------
echo "Step 3: Standardizing chromosome names"
sed 's/^X\t/chrX\t/' "${NOHEADER}" > "${STD_CHR}"
echo "  Records: $(wc -l < "${STD_CHR}")"

# -----------------------------------------------------------------------------
# Step 4: Intersect motif coordinates with cytoband intervals
# -wa retains original motif coordinates
# -wb appends overlapping cytoband record
# -----------------------------------------------------------------------------
echo "Step 4: Intersecting motifs with cytobands"
bedtools intersect -a "${STD_CHR}" -b "${CYTO_TXT}" -wa -wb > "${INTERSECTED}"
echo "  Records: $(wc -l < "${INTERSECTED}")"

# -----------------------------------------------------------------------------
# Step 5: Format final output table
# Input columns from intersection:
#   1 chromosome (motif), 2 start, 3 end, 4 motif,
#   5 chromosome (cytoBand), 6 start (cytoBand), 7 end (cytoBand),
#   8 cytoband name
# Output columns: Motif, Chromosome, Start, End, Cytoband
# -----------------------------------------------------------------------------
echo "Step 5: Formatting final table"
printf "Motif\tChromosome\tStart\tEnd\tCytoband\n" > "${FINAL}"
awk -F'\t' '{print $4"\t"$1"\t"$2"\t"$3"\t"$8}' "${INTERSECTED}" >> "${FINAL}"
echo "  Final table: ${FINAL}"
echo "  Total records (including header): $(wc -l < "${FINAL}")"

# -----------------------------------------------------------------------------
# Step 6: Summary
# -----------------------------------------------------------------------------
echo "Step 6: Summary of results"
echo "  Motifs present:"
awk -F'\t' 'NR>1 {print $1}' "${FINAL}" | sort -u | sed 's/^/    /'
echo "  Per-motif counts:"
awk -F'\t' 'NR>1 {c[$1]++} END {for (m in c) printf "    %-6s %d\n", m, c[m]}' "${FINAL}" | sort
echo "  Chromosome distribution:"
awk -F'\t' 'NR>1 {c[$2]++} END {for (chr in c) printf "    %-8s %d\n", chr, c[chr]}' "${FINAL}" | sort

echo "Done."
