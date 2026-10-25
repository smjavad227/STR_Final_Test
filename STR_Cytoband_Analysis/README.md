# STR Cytoband Annotation Pipeline

Reproducible pipeline for annotating eleven conserved non-CG short tandem
repeat (STR) motifs with chromosomal cytoband locations based on the human
reference genome assembly GRCh38/hg38.

## Overview

The pipeline performs three main tasks:

1. Consolidates motif coordinates from two source files into a single
   four-column BED file containing all eleven motifs.
2. Downloads and prepares the UCSC cytoBand reference for GRCh38.
3. Intersects motif coordinates with cytoband intervals and formats the
   final five-column table for use as a supplementary table in a manuscript.

## Requirements

- Linux or macOS environment (tested on Ubuntu 22.04)
- bedtools v2.30.0 or later
- Python 3.10 or later
- wget and standard command-line tools (grep, awk, sed, sort)

## Installation

Option A: Conda environment

    conda env create -f environment.yml
    conda activate str_cytoband

Option B: Manual installation

    conda install -c bioconda bedtools

## Input Data

Place the two source files in data/input/:

- 11_motifs_coordinates_final.bed
  Four-column, tab-separated file with header.
  Contains coordinates for eight of the eleven motifs.

- Xa_STRs_corrected.tsv
  Tab-separated file containing additional motifs, including
  TATC, AGAT, and AAATT. Only these three motifs are extracted.

## Usage

From the repository root, run:

    bash run.sh

The pipeline performs input preparation and annotation sequentially.
The final output is written to:

    data/output/motifs_11_paper_final_all.tsv

## Output

The final table contains five columns:

    Motif    Chromosome    Start    End    Cytoband

Each row corresponds to a single occurrence of one of the eleven motifs,
annotated with the UCSC cytoband that overlaps its coordinates.

## Reproducibility

All scripts use relative paths resolved from the script location.
The pipeline is independent of the working directory and can be executed
on any system with the required dependencies. Running the pipeline on the
provided input data reproduces the final table deterministically.

## License

MIT