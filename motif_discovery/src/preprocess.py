#!/usr/bin/env python3
"""
preprocess.py
Utilities for reading, cleaning, and preparing biological sequences for motif discovery.

Author: Ramesh Chandra Soren (2022CSB086)
Supervisor: Dr. Surajeet Ghosh, IIEST Shibpur
"""

import re
import os
from collections import Counter


def read_fasta(filepath):
    """
    Parse a FASTA file into a list of (header, sequence) tuples.
    Handles multi-line sequences and IUPAC ambiguity codes.
    """
    sequences = []
    current_header = None
    current_seq = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line.upper())

    if current_header is not None:
        sequences.append((current_header, ''.join(current_seq)))

    return sequences


def clean_sequence(seq, allowed='ACGT', replace_unknown='N'):
    """Remove or replace characters not in the allowed set."""
    cleaned = []
    for ch in seq.upper():
        if ch in allowed:
            cleaned.append(ch)
        else:
            # Skip ambiguous codes like N, R, Y, etc.
            pass
    return ''.join(cleaned)


def filter_sequences(sequences, min_length=50, max_length=5000, allowed='ACGT'):
    """
    Filter and clean sequences:
    - Remove too-short or too-long sequences
    - Remove non-ACGT characters
    - Warn about highly repetitive or low-complexity sequences
    """
    filtered = []
    stats = {'total': len(sequences), 'passed': 0, 'too_short': 0,
             'too_long': 0, 'low_complexity': 0}

    for header, seq in sequences:
        cleaned = clean_sequence(seq, allowed)

        if len(cleaned) < min_length:
            stats['too_short'] += 1
            continue
        if len(cleaned) > max_length:
            cleaned = cleaned[:max_length]  # Truncate to max

        # Low-complexity filter: if one nucleotide > 80%
        counts = Counter(cleaned)
        max_freq = max(counts.values()) / len(cleaned)
        if max_freq > 0.80:
            stats['low_complexity'] += 1
            continue

        filtered.append((header, cleaned))
        stats['passed'] += 1

    return filtered, stats


def compute_background_frequencies(sequences):
    """Compute background nucleotide frequencies from all sequences."""
    total = Counter()
    for _, seq in sequences:
        total.update(seq)
    n = sum(total.values())
    bg = {nuc: total.get(nuc, 0) / n for nuc in 'ACGT'}
    return bg


def write_fasta(sequences, filepath):
    """Write sequences to a FASTA file."""
    with open(filepath, 'w') as f:
        for header, seq in sequences:
            f.write(f'>{header}\n')
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')


def preprocess_dataset(input_path, output_path, min_length=50, max_length=2000):
    """Full preprocessing pipeline for a FASTA file."""
    print(f"\n[Preprocessing] {os.path.basename(input_path)}")
    raw = read_fasta(input_path)
    print(f"  Raw sequences: {len(raw)}")

    filtered, stats = filter_sequences(raw, min_length=min_length, max_length=max_length)
    print(f"  Passed filter: {stats['passed']}")
    print(f"  Too short: {stats['too_short']}, Low complexity: {stats['low_complexity']}")

    bg = compute_background_frequencies(filtered)
    print(f"  Background freqs: A={bg['A']:.3f} C={bg['C']:.3f} "
          f"G={bg['G']:.3f} T={bg['T']:.3f}")

    write_fasta(filtered, output_path)
    print(f"  Saved to: {output_path}")

    return filtered, bg, stats


if __name__ == "__main__":
    import glob
    os.makedirs("../data/processed", exist_ok=True)

    datasets = glob.glob("../data/raw/*.fasta")
    for path in sorted(datasets):
        name = os.path.basename(path)
        out = f"../data/processed/{name}"
        try:
            seqs, bg, stats = preprocess_dataset(path, out)
        except Exception as e:
            print(f"  [WARNING] Could not process {name}: {e}")

    print("\nPreprocessing complete.")
