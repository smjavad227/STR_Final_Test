#!/usr/bin/env bash
# =============================================================================
# 01_prepare_input.sh
#
# Consolidates coordinates for all eleven conserved non-CG STR motifs into
# a single four-column, tab-separated file with a header.
#
# Sources:
#   data/input/11_motifs_coordinates_final.bed
#     Contains coordinates for multiple motifs. Only the eight target motifs
#     are extracted: ACAC, TTTAA, TCTT, GAGG, AAGA, ATCT, ATAG, TATG.
#
#   data/input/Xa_STRs_corrected.tsv
#     Contains additional motif records including TATC, AGAT, and AAATT.
#     Only these three motifs are extracted from this file.
#
# Output:
#   data/intermediate/all_11_motifs_coordinates.bed
#     Four-column file with header: #chr, start, end, motif.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

SRC_8="data/input/11_motifs_coordinates_final.bed"
SRC_3="data/input/Xa_STRs_corrected.tsv"
OUT="data/intermediate/all_11_motifs_coordinates.bed"
TMP_8="data/intermediate/tmp_8motifs.bed"
TMP_3="data/intermediate/tmp_3motifs.bed"

mkdir -p data/intermediate

# Validate source files
if [[ ! -s "${SRC_8}" ]]; then
    echo "ERROR: missing source file: ${SRC_8}"
    exit 1
fi
if [[ ! -s "${SRC_3}" ]]; then
    echo "ERROR: missing source file: ${SRC_3}"
    exit 1
fi

echo "Step 1a: Extracting eight target motifs from ${SRC_8}"
# The source file may contain additional motifs beyond the eight target ones.
# Filter explicitly for the eight target motifs to guarantee the final set.
# Skip the header line first, then filter column 4 by exact match.
tail -n +2 "${SRC_8}" \
    | awk -F'\t' 'BEGIN {OFS="\t"}
        $4=="ACAC" || $4=="TTTAA" || $4=="TCTT" || $4=="GAGG" \
        || $4=="AAGA" || $4=="ATCT" || $4=="ATAG" || $4=="TATG"' \
    > "${TMP_8}"
echo "  Records: $(wc -l < "${TMP_8}")"

echo "Step 1b: Extracting three additional motifs from ${SRC_3}"
# Column 9 in the source file corresponds to the motif sequence.
# Convert to uppercase for consistency with the other source file.
grep -i -w -E "tatc|agat|aaatt" "${SRC_3}" \
    | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"toupper($9)}' > "${TMP_3}"
echo "  Records: $(wc -l < "${TMP_3}")"

echo "Step 1c: Merging into a single file"
printf "#chr\tstart\tend\tmotif\n" > "${OUT}"
cat "${TMP_8}" "${TMP_3}" | sort -u >> "${OUT}"

echo "  Output: ${OUT}"
echo "  Total records (including header): $(wc -l < "${OUT}")"
echo "  Motifs present:"
awk -F'\t' 'NR>1 {print $4}' "${OUT}" | sort -u | sed 's/^/    /'
