#!/usr/bin/env bash
# =============================================================================
# run.sh
#
# Entry point for the STR cytoband annotation pipeline.
# Executes input preparation followed by cytoband annotation.
#
# Usage:
#   bash run.sh
# =============================================================================

set -euo pipefail

# Resolve repository root from the script location (path-independent)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# Ensure output directories exist
mkdir -p data/intermediate data/output data/reference results

# Step 1: Prepare consolidated input file with all eleven motifs
bash scripts/01_prepare_input.sh

# Step 2: Annotate motifs with cytoband locations
bash scripts/02_annotate_motifs.sh

echo "Pipeline complete."
