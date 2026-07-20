import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from config.paths import *

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
evolutionary_conservation_analysis.py

Phase 4: Cross‑species evolutionary conservation analysis of STR motifs.

This script:
    - Loads combined STR data from five primate species (human, chimpanzee,
      gorilla, orangutan, macaque).
    - For each motif, calculates conservation score based on frequency
      consistency across species.
    - Identifies highly conserved non‑CG motifs.
    - Outputs tables and figures for Supplementary Table S3 and Figure 3.

Input:
    all_primate_strs_combined.tsv – combined STR table from all species
        (must contain columns: 'species' and 'consensus_sequence')

Outputs:
    conserved_motifs.csv            – all motifs with conservation scores
    forbidden_motifs_candidates.csv – top conserved non‑CG candidates
    conservation_scores_distribution.png
    forbidden_motifs_heatmap.png
    motif_length_vs_conservation.png
    evolutionary_conservation_summary.txt

Written for:
"Distinct evolutionary signatures shape depletion and conservation of
short tandem repeat motifs in the human genome"
"""

import os
import sys
import argparse
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# =============================================================================
#  Core analysis functions
# =============================================================================

def load_data(file_path):
    """Load the combined multi‑species STR table."""
    df = pd.read_csv(file_path, sep='\t')
    return df


def compute_conservation_scores(df, min_species=2, min_total_count=10):
    """
    Compute conservation score for each motif based on frequency variance
    across species.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain columns 'species' and 'consensus_sequence'.
    min_species : int
        Minimum number of species in which a motif must appear.
    min_total_count : int
        Minimum total occurrence count to be considered.

    Returns
    -------
    list of dict
        Each dict contains: motif, length, species_count, mean_frequency,
        std_frequency, conservation_score, coefficient_of_variation.
    """
    # Count motifs per species
    species_counter = {}
    for species in df['species'].unique():
        sub = df[df['species'] == species]
        species_counter[species] = Counter(sub['consensus_sequence'].dropna())

    all_motifs = set()
    for cnt in species_counter.values():
        all_motifs.update(cnt.keys())

    results = []

    for motif in all_motifs:
        present_species = []
        freq_values = []

        for sp, cnt in species_counter.items():
            count = cnt.get(motif, 0)
            if count > 0:
                present_species.append(sp)
                total_motifs = sum(cnt.values())
                freq = count / total_motifs if total_motifs > 0 else 0.0
                freq_values.append(freq)

        total_count = sum(cnt.get(motif, 0) for cnt in species_counter.values())

        if len(present_species) >= min_species and total_count >= min_total_count:
            mean_freq = np.mean(freq_values)
            std_freq = np.std(freq_values)
            cv = std_freq / mean_freq if mean_freq > 0 else 0.0
            conservation_score = 1.0 - cv  # lower CV → higher conservation
            conservation_score = max(0.0, min(1.0, conservation_score))  # clamp

            results.append({
                'motif': motif,
                'length': len(motif),
                'species_count': len(present_species),
                'mean_frequency': mean_freq,
                'std_frequency': std_freq,
                'conservation_score': conservation_score,
                'coefficient_of_variation': cv,
                'species': ', '.join(sorted(present_species))
            })

    # Sort by conservation score descending
    results.sort(key=lambda x: x['conservation_score'], reverse=True)
    return results


def identify_forbidden_candidates(conserved_list, min_species=3, min_score=0.7):
    """
    Filter for conserved non‑CG motifs that meet strict criteria.
    (The actual CG‑filtering is applied externally; here we just keep those
    with high conservation.)
    """
    candidates = [x for x in conserved_list
                  if x['species_count'] >= min_species
                  and x['conservation_score'] >= min_score]
    return candidates


# =============================================================================
#  Plotting functions
# =============================================================================

def plot_conservation_distribution(conserved_list, out_dir):
    """Histogram of conservation scores."""
    scores = [x['conservation_score'] for x in conserved_list]
    plt.figure(figsize=(10, 6))
    plt.hist(scores, bins=30, alpha=0.7, color='skyblue', edgecolor='black')
    plt.xlabel('Conservation score')
    plt.ylabel('Number of motifs')
    plt.title('Distribution of conservation scores across all motifs')
    plt.grid(alpha=0.3)
    plt.savefig(os.path.join(out_dir, 'conservation_scores_distribution.png'),
                dpi=300, bbox_inches='tight')
    plt.close()


def plot_heatmap(df, candidates, out_dir):
    """
    Heatmap of motif frequencies across species for the top conserved motifs.
    """
    if not candidates:
        return
    # Take top 10 candidates
    top_motifs = [x['motif'] for x in candidates[:10]]
    species_list = sorted(df['species'].unique())

    freq_matrix = np.zeros((len(top_motifs), len(species_list)))
    for i, motif in enumerate(top_motifs):
        for j, sp in enumerate(species_list):
            sub = df[df['species'] == sp]
            total = len(sub)
            if total == 0:
                continue
            cnt = Counter(sub['consensus_sequence'].dropna())
            freq_matrix[i, j] = cnt.get(motif, 0) / total

    plt.figure(figsize=(12, 8))
    sns.heatmap(freq_matrix, annot=True, fmt='.3f',
                xticklabels=species_list, yticklabels=top_motifs,
                cmap='YlOrRd', cbar_kws={'label': 'Frequency'})
    plt.xlabel('Species')
    plt.ylabel('Motif')
    plt.title('Frequency of top conserved motifs across species')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'forbidden_motifs_heatmap.png'),
                dpi=300, bbox_inches='tight')
    plt.close()


def plot_length_vs_conservation(conserved_list, out_dir):
    """Scatter plot of motif length vs conservation score."""
    lengths = [x['length'] for x in conserved_list]
    scores = [x['conservation_score'] for x in conserved_list]

    plt.figure(figsize=(10, 6))
    plt.scatter(lengths, scores, alpha=0.5, color='green')
    plt.xlabel('Motif length (bp)')
    plt.ylabel('Conservation score')
    plt.title('Motif length vs conservation score')
    plt.grid(alpha=0.3)

    # Linear regression line
    if len(lengths) > 1:
        z = np.polyfit(lengths, scores, 1)
        p = np.poly1d(z)
        x_sorted = sorted(lengths)
        plt.plot(x_sorted, p(x_sorted), 'r--', alpha=0.8, label=f'slope={z[0]:.3f}')
        plt.legend()

    plt.savefig(os.path.join(out_dir, 'motif_length_vs_conservation.png'),
                dpi=300, bbox_inches='tight')
    plt.close()


# =============================================================================
#  Command‑line interface
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyse evolutionary conservation of STR motifs."
    )
    parser.add_argument('--input', required=True,
                        help='Path to all_primate_strs_combined.tsv')
    parser.add_argument('--output_dir', required=True,
                        help='Directory for output tables and figures')
    parser.add_argument('--log', default=None,
                        help='Optional log file')
    return parser.parse_args()


def main():
    args = parse_args()

    if args.log:
        sys.stdout = open(args.log, 'w', encoding='utf-8')

    print("=" * 70)
    print("EVOLUTIONARY CONSERVATION ANALYSIS")
    print("=" * 70)

    # 1. Load data
    print(f"\n[1] Loading: {args.input}")
    df = load_data(args.input)
    print(f"    Total STRs: {len(df):,}")
    print(f"    Species:    {df['species'].unique().tolist()}")

    # 2. Compute conservation scores
    print("\n[2] Computing conservation scores...")
    conserved = compute_conservation_scores(df, min_species=2, min_total_count=10)
    print(f"    Conserved motifs found: {len(conserved):,}")

    # 3. Identify candidates
    candidates = identify_forbidden_candidates(conserved, min_species=3, min_score=0.7)
    print(f"    High‑confidence candidates: {len(candidates)}")

    # 4. Save tables
    os.makedirs(args.output_dir, exist_ok=True)
    tbl_dir = os.path.join(args.output_dir, 'tables')
    os.makedirs(tbl_dir, exist_ok=True)

    cons_df = pd.DataFrame(conserved)
    cons_df.to_csv(os.path.join(tbl_dir, 'conserved_motifs.csv'), index=False)
    if candidates:
        cand_df = pd.DataFrame(candidates)
        cand_df.to_csv(os.path.join(tbl_dir, 'forbidden_motifs_candidates.csv'), index=False)

    print("\n[4] Saved tables to:", tbl_dir)

    # 5. Generate figures
    fig_dir = os.path.join(args.output_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)

    plot_conservation_distribution(conserved, fig_dir)
    plot_heatmap(df, candidates, fig_dir)
    plot_length_vs_conservation(conserved, fig_dir)
    print("    Figures saved to:", fig_dir)

    # 6. Summary report
    report = f"""
================================================================
EVOLUTIONARY CONSERVATION ANALYSIS – SUMMARY
================================================================
Input:           {args.input}
Total STRs:      {len(df):,}
Species:         {', '.join(df['species'].unique())}

Conserved motifs (min 2 species, ≥10 occurrences): {len(conserved):,}
Average conservation score: {np.mean([x['conservation_score'] for x in conserved]):.3f}
Range:           {min([x['conservation_score'] for x in conserved]):.3f} – {max([x['conservation_score'] for x in conserved]):.3f}

High‑confidence candidates (≥3 species, score ≥0.7): {len(candidates)}

Top 10 candidates:
"""
    for i, item in enumerate(candidates[:10], 1):
        report += f"  {i:2d}. {item['motif']}  (len={item['length']})  species={item['species_count']}  score={item['conservation_score']:.3f}\n"

    report += f"""
Output:
  Tables:  {tbl_dir}
  Figures: {fig_dir}
================================================================
"""
    with open(
        os.path.join(args.output_dir, "evolutionary_conservation_summary.txt"),
        "w",
        encoding="utf-8"
    ) as f:
        f.write(report)

    print("\n" + report)
    print("Evolutionary conservation analysis completed.")


if __name__ == "__main__":
    main()