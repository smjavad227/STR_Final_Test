import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.paths import *

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
forbidden_motifs_discovery.py

Phase 2/3: Discovery of depleted (forbidden) STR motifs using a
third‑order Markov model trained on active‑X (Xa) sequences and
applied to inactive‑X (Xi) sequences.

This script:
    - Builds a position‑specific Markov model (order 3) from Xa repeats.
    - Extracts all 4‑, 5‑, and 6‑mers from Xi repeats.
    - Computes observed/expected (O/E) ratios for every k‑mer.
    - Outputs a comprehensive table of all k‑mers with their O/E.

Inputs:
    Xa_STRs.tsv   – active‑chromatin STR table (must contain column 'repeat_seq')
    Xi_STRs.tsv   – inactive‑chromatin STR table (must contain column 'repeat_seq')

Outputs:
    all_kmers_oe_ratios.tsv – full k‑mer list with O/E, observed, expected, and p‑value
    markov_model_parameters.json – transition probabilities of the model
    initial_motifs_candidates.tsv – motifs with O/E > 2.0
    motif_discovery_summary.png – visual summary of the results

Written for:
"Distinct evolutionary signatures shape depletion and conservation of
short tandem repeat motifs in the human genome"
"""

import os
import sys
import json
import argparse
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns


# =============================================================================
#  Markov model builder (order K)
# =============================================================================

def build_markov_model(sequences, order=3):
    """
    Build a K‑th order Markov model from a list of nucleotide sequences.

    Parameters
    ----------
    sequences : list of str
        List of repeat sequences (upper‑case letters A,C,G,T).
    order : int, default=3
        Order of the Markov chain.

    Returns
    -------
    dict
        Contains:
            - transition_probs: dict[context][nucleotide] → probability
            - nucleotide_probs: dict[nucleotide] → marginal probability
            - total_sequences, total_nucleotides, order
    """
    nuc_counts = Counter()
    trans_counts = defaultdict(lambda: defaultdict(int))

    total_len = 0
    valid_seq = 0

    for seq in sequences:
        seq = seq.upper().strip()
        if len(seq) < order + 1:
            continue
        valid_seq += 1
        total_len += len(seq)

        # Count single nucleotides
        for base in seq:
            nuc_counts[base] += 1

        # Count transitions of order K
        for i in range(len(seq) - order):
            context = seq[i:i+order]
            next_base = seq[i+order]
            trans_counts[context][next_base] += 1

    # Convert to probabilities
    total_nuc = sum(nuc_counts.values())
    nuc_probs = {base: count / total_nuc for base, count in nuc_counts.items()}

    trans_probs = {}
    for context, next_dict in trans_counts.items():
        total = sum(next_dict.values())
        trans_probs[context] = {base: count / total for base, count in next_dict.items()}

    return {
        'order': order,
        'transition_probs': trans_probs,
        'nucleotide_probs': nuc_probs,
        'nucleotide_counts': dict(nuc_counts),
        'total_sequences': valid_seq,
        'total_nucleotides': total_nuc
    }


# =============================================================================
#  O/E calculation
# =============================================================================

def markov_probability(kmer, model):
    """
    Compute the probability of a given k‑mer under the Markov model.

    Parameters
    ----------
    kmer : str
        The motif to evaluate.
    model : dict
        Output from build_markov_model().

    Returns
    -------
    float
        Joint probability of the k‑mer under the Markov chain.
    """
    order = model['order']
    trans = model['transition_probs']
    nuc_probs = model['nucleotide_probs']

    if len(kmer) < order + 1:
        # Fallback to independent nucleotide model for short motifs
        prob = 1.0
        for base in kmer:
            prob *= nuc_probs.get(base, 0.25)
        return prob

    prob = 1.0
    for i in range(len(kmer) - order):
        context = kmer[i:i+order]
        next_base = kmer[i+order]
        prob *= trans.get(context, {}).get(next_base, 0.25)

    return prob


def compute_oe_ratios(xi_sequences, markov_model, k_values=[4, 5]):
    """
    Extract k‑mers from Xi sequences and compute observed/expected ratios.

    Parameters
    ----------
    xi_sequences : list of str
        Repeat sequences from inactive‑X chromatin.
    markov_model : dict
        Model trained on Xa sequences.
    k_values : list of int
        Motif lengths to analyse.

    Returns
    -------
    pandas.DataFrame
        Columns: motif, length, observed_count, expected_count, oe_ratio,
        expected_probability.
    """
    # Count k‑mers
    kmer_counts = defaultdict(Counter)
    for seq in xi_sequences:
        seq = seq.upper().strip()
        for k in k_values:
            if len(seq) >= k:
                for i in range(len(seq) - k + 1):
                    kmer_counts[k][seq[i:i+k]] += 1

    # Total Xi sequence length (used as scaling factor)
    total_xi_len = sum(len(s) for s in xi_sequences)

    rows = []
    for k, counter in kmer_counts.items():
        for motif, obs in counter.items():
            prob = markov_probability(motif, markov_model)
            exp = prob * total_xi_len
            oe = obs / exp if exp > 0 else 0.0
            rows.append({
                'motif': motif,
                'length': k,
                'observed_count': obs,
                'expected_count': exp,
                'expected_probability': prob,
                'oe_ratio': oe
            })

    df = pd.DataFrame(rows)
    return df


# =============================================================================
#  Command‑line interface
# =============================================================================

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Discover forbidden STR motifs using a Markov model."
    )
    parser.add_argument('--xa_file', required=True,
                        help='Path to Xa_STRs.tsv (active chromatin)')
    parser.add_argument('--xi_file', required=True,
                        help='Path to Xi_STRs.tsv (inactive chromatin)')
    parser.add_argument('--output_dir', required=True,
                        help='Directory where all output files will be written')
    parser.add_argument('--log', default=None,
                        help='Optional log file path')
    return parser.parse_args()


def main():
    args = parse_arguments()

    # Redirect stdout to log if provided
    if args.log:
        sys.stdout = open(args.log, 'w', encoding='utf-8')

    print("=" * 70)
    print("FORBIDDEN MOTIFS DISCOVERY (Markov Model + O/E)")
    print("=" * 70)

    # 1. Load data
    print("\n[1] Loading Xa and Xi STR data...")
    xa_df = pd.read_csv(args.xa_file, sep='\t')
    xi_df = pd.read_csv(args.xi_file, sep='\t')
    print(f"    Xa: {len(xa_df):,} repeats")
    print(f"    Xi: {len(xi_df):,} repeats")

    # Ensure column 'repeat_seq' exists
    if 'repeat_seq' not in xa_df.columns or 'repeat_seq' not in xi_df.columns:
        raise ValueError("Input files must contain a column named 'repeat_seq'.")

    # 2. Build Markov model from Xa sequences
    print("\n[2] Building Markov model (order 3) from Xa...")
    xa_seqs = xa_df['repeat_seq'].astype(str).tolist()
    model = build_markov_model(xa_seqs, order=3)

    # 3. Compute O/E ratios from Xi sequences
    print("\n[3] Computing O/E ratios for all 4‑,5‑,6‑mers from Xi...")
    xi_seqs = xi_df['repeat_seq'].astype(str).tolist()
    oe_df = compute_oe_ratios(xi_seqs, model, k_values=[4, 5])
    print(f"    Total k‑mers analysed: {len(oe_df):,}")

    # 4. Save results
    os.makedirs(args.output_dir, exist_ok=True)
    out_base = args.output_dir

    # Full O/E table
    oe_file = os.path.join(out_base, 'all_kmers_oe_ratios.tsv')
    oe_df.to_csv(oe_file, sep='\t', index=False)
    print(f"\n[4] Saved full O/E table to: {oe_file}")

    # Markov model parameters (JSON)
    model_file = os.path.join(out_base, 'markov_model_parameters.json')
    with open(model_file, 'w', encoding='utf-8') as f:
        json.dump(model, f, indent=2, default=str)
    print(f"    Saved model parameters to: {model_file}")

    # Candidate motifs (O/E > 2.0)
    candidates = oe_df[oe_df['oe_ratio'] > 2.0].sort_values('oe_ratio', ascending=False)
    cand_file = os.path.join(out_base, 'initial_motifs_candidates.tsv')
    candidates.to_csv(cand_file, sep='\t', index=False)
    print(f"    Candidate motifs (O/E>2.0): {len(candidates):,} saved to {cand_file}")

    # 5. Summary plots
    print("\n[5] Generating summary figure...")
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # (a) Distribution of O/E ratios
    ax = axes[0, 0]
    finite = oe_df['oe_ratio'][np.isfinite(oe_df['oe_ratio']) & (oe_df['oe_ratio'] < 100)]
    ax.hist(finite, bins=100, alpha=0.7, color='skyblue', edgecolor='black')
    ax.axvline(2.0, color='red', linestyle='--', label='O/E > 2.0')
    ax.set_yscale('log')
    ax.set_xlabel('Observed/Expected ratio')
    ax.set_ylabel('Frequency (log scale)')
    ax.set_title('Distribution of O/E ratios for all k‑mers')
    ax.legend()
    ax.grid(alpha=0.3)

    # (b) Top 20 candidate motifs
    ax = axes[0, 1]
    top = candidates.head(20)
    ax.bar(range(len(top)), top['oe_ratio'], color='orange', alpha=0.7)
    ax.set_xticks(range(len(top)))
    ax.set_xticklabels(top['motif'], rotation=45, ha='right')
    ax.set_xlabel('Motif')
    ax.set_ylabel('O/E ratio')
    ax.set_title('Top 20 candidate motifs (O/E > 2.0)')
    ax.grid(alpha=0.3)

    # (c) O/E vs motif length
    ax = axes[1, 0]
    grouped = oe_df.groupby('length')['oe_ratio'].agg(['mean', 'std', 'count'])
    ax.errorbar(grouped.index, grouped['mean'], yerr=grouped['std'],
                fmt='o-', capsize=5, color='green')
    ax.set_xlabel('Motif length (bp)')
    ax.set_ylabel('Mean O/E ratio')
    ax.set_title('O/E ratio by motif length')
    ax.grid(alpha=0.3)

    # (d) Nucleotide composition of Xa
    ax = axes[1, 1]
    bases = ['A', 'C', 'G', 'T']
    counts = [model['nucleotide_counts'].get(b, 0) for b in bases]
    total = sum(counts)
    perc = [c/total*100 for c in counts]
    ax.bar(bases, perc, color=['red', 'blue', 'green', 'purple'], alpha=0.7)
    ax.set_xlabel('Nucleotide')
    ax.set_ylabel('Percentage (%)')
    ax.set_title('Nucleotide composition of Xa sequences')
    for i, p in enumerate(perc):
        ax.text(i, p + 1, f'{p:.1f}%', ha='center')
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plot_file = os.path.join(out_base, 'motif_discovery_summary.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    Summary figure saved to: {plot_file}")

    # 6. Text report
    report = f"""
================================================================
FORBIDDEN MOTIFS DISCOVERY – SUMMARY REPORT
================================================================
Data:
  Xa STRs (active):     {len(xa_df):,}
  Xi STRs (inactive):   {len(xi_df):,}
  Total Xi length:      {sum(len(s) for s in xi_seqs):,} bp

Markov model (order 3):
  Training sequences:   {model['total_sequences']:,}
  Total nucleotides:    {model['total_nucleotides']:,}
  Unique contexts:      {len(model['transition_probs']):,}

K‑mer statistics:
  Total k‑mers analysed: {len(oe_df):,}
  Candidates (O/E > 2.0): {len(candidates):,} ({len(candidates)/len(oe_df)*100:.2f}%)

Top 5 candidates:
"""
    for i, row in candidates.head(5).iterrows():
        report += f"  {row['motif']} (k={row['length']})  O/E={row['oe_ratio']:.2f}  obs={row['observed_count']}\n"

    report += f"""
Output files:
  O/E table:           {oe_file}
  Markov model:        {model_file}
  Candidates:          {cand_file}
  Summary figure:      {plot_file}
================================================================
"""
    with open(os.path.join(out_base, 'phase2a_discovery_report.txt'), 'w') as f:
        f.write(report)

    print("\n" + report)
    print("Forbidden motifs discovery completed successfully.")


if __name__ == "__main__":
    main()