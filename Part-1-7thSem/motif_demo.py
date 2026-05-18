#!/usr/bin/env python3
"""
motif_demo.py
- Reads /home/ramesh/Desktop/Major_p/motifs.fasta
- Runs Greedy Motif Search for k in [5,6,7,8]
- Picks the best motifs (lowest Hamming score), writes results to results.txt
- Saves PWM heatmap (pwm.png) and a sequence-logo (logo.png) if possible
"""
import os
import re
from collections import defaultdict
import json

DATA_PATH = os.path.join(os.path.dirname(__file__), 'motifs.fasta')
OUT_PREFIX = os.path.join(os.path.dirname(__file__), 'motif_demo')

# --- Utilities ---

def read_simple_fasta(path):
    seqs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # accept only lines that look like DNA
            if re.fullmatch(r'[ACGTacgt]+', line):
                seqs.append(line.upper())
    return seqs

# Greedy motif search (from notebook)

def build_profile(motifs):
    k = len(motifs[0])
    num_motifs = len(motifs)
    profile = { 'A': [1] * k, 'C': [1] * k, 'G': [1] * k, 'T': [1] * k }
    for dna_string in motifs:
        for i, nucleotide in enumerate(dna_string):
            if nucleotide in profile:
                profile[nucleotide][i] += 1
    total = num_motifs + 4
    for nucleotide in profile:
        for i in range(k):
            profile[nucleotide][i] /= total
    return profile


def score_kmer(kmer, profile):
    score = 1.0
    for i, nucleotide in enumerate(kmer):
        if nucleotide in profile:
            score *= profile[nucleotide][i]
        else:
            return 0.0
    return score


def find_best_kmer(sequence, k, profile):
    best_kmer = sequence[0:k]
    best_score = score_kmer(best_kmer, profile)
    for i in range(1, len(sequence) - k + 1):
        current_kmer = sequence[i:i+k]
        current_score = score_kmer(current_kmer, profile)
        if current_score > best_score:
            best_score = current_score
            best_kmer = current_kmer
    return best_kmer


def get_consensus(profile):
    k = len(profile['A'])
    consensus = ''
    for i in range(k):
        max_prob = -1
        best_nuc = ''
        for nucleotide in ['A','C','G','T']:
            if profile[nucleotide][i] > max_prob:
                max_prob = profile[nucleotide][i]
                best_nuc = nucleotide
        consensus += best_nuc
    return consensus


def calculate_score(motifs, consensus):
    score = 0
    for motif in motifs:
        for i in range(len(motif)):
            if motif[i] != consensus[i]:
                score += 1
    return score


def greedy_motif_search(dna_list, k):
    if not dna_list:
        return [], 0, ''
    best_motifs = [seq[0:k] for seq in dna_list]
    profile = build_profile(best_motifs)
    consensus = get_consensus(profile)
    best_score = calculate_score(best_motifs, consensus)
    first_sequence = dna_list[0]
    for i in range(len(first_sequence) - k + 1):
        seed_kmer = first_sequence[i:i+k]
        current_motifs = [seed_kmer]
        for j in range(1, len(dna_list)):
            profile = build_profile(current_motifs)
            next_motif = find_best_kmer(dna_list[j], k, profile)
            current_motifs.append(next_motif)
        final_profile = build_profile(current_motifs)
        final_consensus = get_consensus(final_profile)
        current_score = calculate_score(current_motifs, final_consensus)
        if current_score < best_score:
            best_score = current_score
            best_motifs = current_motifs
    final_profile = build_profile(best_motifs)
    final_consensus = get_consensus(final_profile)
    return best_motifs, best_score, final_consensus

# PWM builder that returns counts and probabilities

def build_pwm_from_motifs(motifs):
    k = len(motifs[0])
    counts = { 'A':[0]*k, 'C':[0]*k, 'G':[0]*k, 'T':[0]*k }
    for m in motifs:
        for i,ch in enumerate(m):
            counts[ch][i] += 1
    # convert to probabilities
    probs = {n: [ (counts[n][i] / len(motifs)) for i in range(k) ] for n in counts}
    return counts, probs

# Plotting (matplotlib); try logomaker for logo

def save_pwm_heatmap(probs, outpath):
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as e:
        print('matplotlib required for plots. error:', e)
        return
    order = ['A','C','G','T']
    arr = [[probs[n][i] for i in range(len(probs['A']))] for n in order]
    arr = np.array(arr)
    fig, ax = plt.subplots(figsize=(max(4, arr.shape[1]*0.6), 3))
    im = ax.imshow(arr, aspect='auto', cmap='Blues', vmin=0, vmax=1)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order)
    ax.set_xticks(range(arr.shape[1]))
    ax.set_xticklabels([str(i+1) for i in range(arr.shape[1])])
    ax.set_xlabel('Position')
    ax.set_title('PWM (probabilities)')
    fig.colorbar(im, ax=ax, orientation='vertical', fraction=0.02)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def save_sequence_logo(counts, outpath):
    # try logomaker first
    try:
        import pandas as pd
        import numpy as np
        import matplotlib.pyplot as plt
        import logomaker
        df = pd.DataFrame({n: counts[n] for n in ['A','C','G','T']})
        df = df / df.sum(axis=1).max()  # normalize so heights look OK
        plt.figure(figsize=(max(6, df.shape[0]*0.6), 3))
        logo = logomaker.Logo(df)
        plt.title('Sequence logo (approx)')
        plt.xlabel('Position')
        plt.ylabel('Relative freq')
        plt.tight_layout()
        plt.savefig(outpath, dpi=150)
        plt.close()
        return True
    except Exception as e:
        # fallback to stacked bar chart
        try:
            import matplotlib.pyplot as plt
            import numpy as np
        except Exception:
            print('matplotlib not available for logo fallback')
            return False
        order = ['A','C','G','T']
        k = len(counts['A'])
        ind = range(k)
        bottom = [0]*k
        colors = {'A':'#1f77b4','C':'#ff7f0e','G':'#2ca02c','T':'#d62728'}
        fig, ax = plt.subplots(figsize=(max(6, k*0.6), 3))
        for n in order:
            ax.bar(ind, counts[n], bottom=bottom, color=colors[n], label=n)
            bottom = [bottom[i] + counts[n][i] for i in range(k)]
        ax.set_xticks(ind)
        ax.set_xticklabels([str(i+1) for i in ind])
        ax.set_xlabel('Position')
        ax.set_ylabel('Counts')
        ax.set_title('Stacked counts (logo fallback)')
        ax.legend()
        plt.tight_layout()
        plt.savefig(outpath, dpi=150)
        plt.close()
        return True

# --- Main ---

def main():
    seqs = read_simple_fasta(DATA_PATH)
    if not seqs:
        print('No sequences found in', DATA_PATH)
        return
    # try multiple k values and pick best
    results = {}
    best_overall = None
    for k in range(5,9):
        # skip if sequences shorter than k
        if any(len(s) < k for s in seqs):
            continue
        motifs, score, consensus = greedy_motif_search(seqs, k)
        counts, probs = build_pwm_from_motifs(motifs)
        results[k] = { 'motifs': motifs, 'score': score, 'consensus': consensus, 'counts': counts, 'probs': probs }
        if best_overall is None or score < best_overall['score']:
            best_overall = { 'k': k, 'motifs': motifs, 'score': score, 'consensus': consensus, 'counts': counts, 'probs': probs }

    if best_overall is None:
        print('No valid k found for provided sequences')
        return

    out_json = OUT_PREFIX + '_results.json'
    with open(out_json, 'w') as f:
        json.dump(best_overall, f, indent=2)

    # write human-readable results
    out_txt = OUT_PREFIX + '_results.txt'
    with open(out_txt, 'w') as f:
        f.write(f"Best k: {best_overall['k']}\n")
        f.write(f"Score (total mismatches): {best_overall['score']}\n")
        f.write(f"Consensus: {best_overall['consensus']}\n\nMotifs (one per sequence):\n")
        for m in best_overall['motifs']:
            f.write(m + '\n')
    print('Wrote', out_json, 'and', out_txt)

    # save PWM heatmap and logo
    pwm_png = OUT_PREFIX + '_pwm.png'
    save_pwm_heatmap(best_overall['probs'], pwm_png)
    logo_png = OUT_PREFIX + '_logo.png'
    ok = save_sequence_logo(best_overall['counts'], logo_png)
    if ok:
        print('Saved images:', pwm_png, logo_png)
    else:
        print('Saved PWM; logo not created')

if __name__ == '__main__':
    main()
