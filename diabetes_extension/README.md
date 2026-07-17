# STR Diabetes Extension Project

Evolutionary conserved non-CG motifs and their epigenetic dysregulation in type 1 and type 2 diabetes

## Overview

This repository contains the extended analysis of conserved non-CG motifs on chromosome X in pancreatic beta cells, focusing on type 1 diabetes (T1D) and type 2 diabetes (T2D).

Key findings:
- In healthy beta cells, genes near conserved motifs show significantly higher expression (p = 0.00118, FC = 1.76; Bonferroni-adjusted p = 0.00354).
- This advantage is lost in T1D (FC = 1.11, ns).
- In T2D, the pattern reverses (FC = 0.84, ns).
- Four motifs show significant hypomethylation in T2D (p = 0.0034, FC = 0.53).

## Repository Structure

- data/ : Input files
- scripts/ : Analysis scripts (6 Python files)
- results/ : Final figures and tables

## Requirements

Python 3.10+ with: pandas, numpy, scipy, matplotlib, seaborn, intervaltree, python-docx, requests

## Reproduce

pip install -r requirements.txt
cd scripts
python run_27k_from_full_table.py
python extract_t1d_ifna_methylation.py
python final_analysis_complete.py
python final_figures.py
python generate_word_with_full_legends.py

## Results Summary

Supplementary Table S7: Expression
Normal: Target=0.950, Control=0.539, FC=1.76, p=0.00118, Sig=**
T1D: Target=0.266, Control=0.239, FC=1.11, p=1.000, Sig=ns
T2D: Target=0.179, Control=0.212, FC=0.84, p=1.000, Sig=ns

Supplementary Table S8: Methylation (T2D only)
AAGA: Control=0.104, Disease=0.055, FC=0.53, p=0.0034, Sig=**
GAGG: Control=0.104, Disease=0.055, FC=0.53, p=0.0034, Sig=**
TATG: Control=0.104, Disease=0.055, FC=0.53, p=0.0034, Sig=**
ATCT: Control=0.104, Disease=0.055, FC=0.53, p=0.0034, Sig=**

## Figures

Figure 7A: Heatmap of methylation fold changes
Figure 7B: Bar plot with significance stars

## Citation

"Distinct evolutionary signatures shape depletion and conservation of short tandem repeat motifs in the human genome" (manuscript under review)

## License

MIT
