#!/usr/bin/env python3
"""
evaluate.py
Performance evaluation metrics for motif discovery algorithms.
Compares predicted motif sites against true planted positions.

Metrics:
  - Sensitivity (Recall)
  - Positive Predictive Value (PPV / Precision)
  - Nucleotide Correlation Coefficient (nCC)
  - Site-level accuracy
  - Motif similarity (PWM correlation)

Author: Ramesh Chandra Soren (2022CSB086)
Supervisor: Dr. Surajeet Ghosh, IIEST Shibpur
"""

import math
import numpy as np
from pwm import NUCLEOTIDES, NUC_INDEX


# ──────────────────────────────────────────────────────────────────────────────
# Site-level evaluation
# ──────────────────────────────────────────────────────────────────────────────

def nucleotide_overlap(pred_start, true_start, width):
    """Number of nucleotide positions that overlap between prediction and truth."""
    pred_end  = pred_start + width
    true_end  = true_start + width
    overlap_start = max(pred_start, true_start)
    overlap_end   = min(pred_end,   true_end)
    return max(0, overlap_end - overlap_start)


def site_metrics(predicted_sites, true_positions, width, overlap_threshold=0.5):
    """
    predicted_sites : list of dicts {'seq_index': int, 'position': int, ...}
    true_positions  : dict {seq_index: true_start_position}  OR list (indexed)
    width           : motif width
    overlap_threshold: fraction of width that must overlap to count as correct

    Returns dict of sensitivity, ppv, f1, nCC.
    """
    if isinstance(true_positions, list):
        true_dict = {i: p for i, p in enumerate(true_positions)}
    else:
        true_dict = true_positions

    tp_nuc = 0   # true  positive nucleotides
    fp_nuc = 0   # false positive nucleotides
    fn_nuc = 0   # false negative nucleotides

    seq_indices = set(true_dict.keys())
    predicted_dict = {}
    for site in predicted_sites:
        idx = site.get('seq_index', site.get('seq_idx', 0))
        predicted_dict[idx] = site['position']

    required_overlap = int(math.ceil(overlap_threshold * width))

    correct = 0
    for idx in seq_indices:
        true_p = true_dict[idx]
        pred_p = predicted_dict.get(idx, None)

        if pred_p is not None:
            overlap = nucleotide_overlap(pred_p, true_p, width)
            tp_nuc += overlap
            fp_nuc += width - overlap
            fn_nuc += width - overlap
            if overlap >= required_overlap:
                correct += 1
        else:
            fn_nuc += width

    n_pred = len(predicted_dict)
    n_true = len(true_dict)

    sensitivity = tp_nuc / (tp_nuc + fn_nuc + 1e-10)
    ppv         = tp_nuc / (tp_nuc + fp_nuc + 1e-10)
    f1          = 2 * sensitivity * ppv / (sensitivity + ppv + 1e-10)

    # nCC (nucleotide correlation coefficient, simplified)
    nCC = (tp_nuc / (width * n_true + 1e-10)
           * tp_nuc / (width * n_pred + 1e-10)) ** 0.5

    site_accuracy = correct / max(1, n_true)

    return {
        'sensitivity':    round(sensitivity, 4),
        'ppv':            round(ppv, 4),
        'f1':             round(f1, 4),
        'nCC':            round(nCC, 4),
        'site_accuracy':  round(site_accuracy, 4),
        'n_correct_sites': correct,
        'n_predicted':    n_pred,
        'n_true':         n_true,
    }


# ──────────────────────────────────────────────────────────────────────────────
# Motif similarity (compare two PWMs)
# ──────────────────────────────────────────────────────────────────────────────

def pwm_pearson_correlation(pwm1, pwm2):
    """
    Pearson correlation between two PWM frequency matrices.
    PWMs must have the same width.
    Returns correlation in [-1, 1].
    """
    assert pwm1.width == pwm2.width, "PWM widths must match"
    fm1 = pwm1.frequency_matrix()
    fm2 = pwm2.frequency_matrix()

    v1 = np.array([fm1[n][p] for p in range(pwm1.width) for n in NUCLEOTIDES])
    v2 = np.array([fm2[n][p] for p in range(pwm2.width) for n in NUCLEOTIDES])

    corr = np.corrcoef(v1, v2)[0, 1]
    return float(corr)


def motif_similarity_score(pwm1, pwm2):
    """
    Motif similarity: max Pearson r over all alignments/shifts.
    Only valid when widths match; otherwise use sliding approach.
    """
    if pwm1.width == pwm2.width:
        return pwm_pearson_correlation(pwm1, pwm2)

    # Slide shorter over longer
    shorter, longer = (pwm1, pwm2) if pwm1.width < pwm2.width else (pwm2, pwm1)
    best = -1.0
    for offset in range(longer.width - shorter.width + 1):
        # Build a sub-PWM from longer at this offset
        from pwm import PWM
        sub = PWM(shorter.width, background=longer.background)
        sub.counts = longer.counts[offset:offset + shorter.width].copy()
        sub.nsites = longer.nsites
        sub.finalize()
        r = pwm_pearson_correlation(shorter, sub)
        best = max(best, r)
    return best


# ──────────────────────────────────────────────────────────────────────────────
# Read true positions from FASTA headers
# ──────────────────────────────────────────────────────────────────────────────

def read_true_positions_from_fasta(fasta_path):
    """
    Parse true planted positions from synthetic FASTA headers.
    Expected header format: >seq_001 planted_motif_pos=42
    Returns dict {seq_index: position}.
    """
    true_pos = {}
    import re
    with open(fasta_path) as f:
        idx = 0
        for line in f:
            if line.startswith('>'):
                m = re.search(r'planted_motif_pos=(\d+)', line)
                if m:
                    true_pos[idx] = int(m.group(1))
                idx += 1
    return true_pos


# ──────────────────────────────────────────────────────────────────────────────
# Pretty print results table
# ──────────────────────────────────────────────────────────────────────────────

def print_evaluation_table(results_list):
    """
    results_list: list of dicts with keys
        'algorithm', 'dataset', 'sensitivity', 'ppv', 'f1', 'nCC',
        'site_accuracy', 'information_content'
    """
    cols = ['algorithm', 'dataset', 'sensitivity', 'ppv', 'f1',
            'nCC', 'site_accuracy', 'information_content']
    widths = [12, 20, 12, 8, 8, 8, 14, 20]
    header = '  '.join(f'{c.upper():<{w}}' for c, w in zip(cols, widths))

    print("\n" + "=" * len(header))
    print(header)
    print("=" * len(header))

    for r in results_list:
        row = '  '.join(
            f'{str(r.get(c, "")):<{w}}' for c, w in zip(cols, widths)
        )
        print(row)
    print("=" * len(header))
