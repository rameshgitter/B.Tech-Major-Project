#!/usr/bin/env python3
"""
generate_synthetic_data.py
Generates synthetic benchmark datasets with planted motifs for algorithm validation.
Used in: Biological Sequence Motif Discovery - B.Tech Major Project
Author: Ramesh Chandra Soren (2022CSB086)
Supervisor: Dr. Surajeet Ghosh, IIEST Shibpur
"""

import random
import os

random.seed(42)

NUCLEOTIDES = ['A', 'C', 'G', 'T']

def random_sequence(length):
    return ''.join(random.choices(NUCLEOTIDES, weights=[0.25, 0.25, 0.25, 0.25], k=length))

def mutate_motif(motif, mutation_rate=0.1):
    """Introduce mutations into a motif instance."""
    result = list(motif)
    for i in range(len(result)):
        if random.random() < mutation_rate:
            result[i] = random.choice([n for n in NUCLEOTIDES if n != result[i]])
    return ''.join(result)

def plant_motif(seq_length, motif, mutation_rate=0.1):
    """Plant a motif instance at a random position in a random background."""
    background = random_sequence(seq_length)
    pos = random.randint(0, seq_length - len(motif))
    mutated = mutate_motif(motif, mutation_rate)
    planted = background[:pos] + mutated + background[pos + len(motif):]
    return planted, pos

def generate_benchmark(filename, motif, n_sequences=20, seq_length=100,
                        mutation_rate=0.15, description=""):
    """Generate a FASTA file with planted motifs."""
    records = []
    true_positions = []
    for i in range(n_sequences):
        seq, pos = plant_motif(seq_length, motif, mutation_rate)
        records.append((f"seq_{i+1:03d}", seq, pos))
        true_positions.append(pos)

    with open(filename, 'w') as f:
        for name, seq, pos in records:
            f.write(f">{name} planted_motif_pos={pos}\n")
            # Write in 80-char lines
            for j in range(0, len(seq), 80):
                f.write(seq[j:j+80] + '\n')

    return true_positions

if __name__ == "__main__":
    os.makedirs("../data/raw", exist_ok=True)

    # Dataset 1: TATA-box like motif (yeast promoters)
    tata_motif = "TATATAAG"
    positions1 = generate_benchmark(
        "../data/raw/synthetic_tata_benchmark.fasta",
        tata_motif, n_sequences=20, seq_length=200, mutation_rate=0.1,
        description="TATA-box motif benchmark"
    )
    print(f"[Dataset 1] TATA-box benchmark: motif='{tata_motif}', {len(positions1)} sequences")

    # Dataset 2: E-box motif (CACGTG) - MYC binding
    ebox_motif = "CACGTG"
    positions2 = generate_benchmark(
        "../data/raw/synthetic_ebox_benchmark.fasta",
        ebox_motif, n_sequences=20, seq_length=150, mutation_rate=0.05,
        description="E-box motif benchmark"
    )
    print(f"[Dataset 2] E-box benchmark: motif='{ebox_motif}', {len(positions2)} sequences")

    # Dataset 3: NF-kB binding site
    nfkb_motif = "GGGACTTTCC"
    positions3 = generate_benchmark(
        "../data/raw/synthetic_nfkb_benchmark.fasta",
        nfkb_motif, n_sequences=25, seq_length=250, mutation_rate=0.2,
        description="NF-kB motif benchmark"
    )
    print(f"[Dataset 3] NF-kB benchmark: motif='{nfkb_motif}', {len(positions3)} sequences")

    # Dataset 4: Zinc finger binding (harder, longer motif)
    znf_motif = "GCGTGGCG"
    positions4 = generate_benchmark(
        "../data/raw/synthetic_znf_benchmark.fasta",
        znf_motif, n_sequences=30, seq_length=300, mutation_rate=0.25,
        description="Zinc finger motif benchmark"
    )
    print(f"[Dataset 4] Zinc finger benchmark: motif='{znf_motif}', {len(positions4)} sequences")

    print("\nAll synthetic datasets generated successfully.")
