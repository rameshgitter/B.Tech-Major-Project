#!/usr/bin/env python3
"""
visualize.py
Publication-quality visualisations for motif discovery results.
  - Sequence logo (information-content weighted)
  - PWM heatmap
  - Gibbs convergence plot
  - Algorithm comparison bar chart

Author: Ramesh Chandra Soren (2022CSB086)
Supervisor: Dr. Surajeet Ghosh, IIEST Shibpur
"""

import os
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.font_manager import FontProperties

COLORS = {'A': '#2ecc71', 'C': '#3498db', 'G': '#f1c40f', 'T': '#e74c3c'}
NUCLEOTIDES = ['A', 'C', 'G', 'T']
NUC_INDEX = {n: i for i, n in enumerate(NUCLEOTIDES)}


# ──────────────────────────────────────────────────────────────────────────────
# 1. Sequence Logo
# ──────────────────────────────────────────────────────────────────────────────

def plot_sequence_logo(pwm, title="Sequence Logo", output_path=None,
                       figsize=(10, 3)):
    """
    Draw a sequence logo from a PWM object.
    Letter height = information content × frequency.
    """
    freq_matrix = pwm.frequency_matrix()
    width = pwm.width

    # Information content per position (bits)
    bg = np.array([pwm.background.get(n, 0.25) for n in NUCLEOTIDES])
    ic_per_pos = []
    for p in range(width):
        freqs = np.array([freq_matrix[n][p] for n in NUCLEOTIDES])
        freqs = np.clip(freqs, 1e-10, 1.0)
        ic = 2.0 + np.sum(freqs * np.log2(freqs))   # max=2 for DNA
        ic_per_pos.append(max(0, ic))

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(0, width)
    ax.set_ylim(0, 2.1)
    ax.set_xlabel("Position", fontsize=12)
    ax.set_ylabel("Information Content (bits)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xticks(np.arange(0.5, width, 1))
    ax.set_xticklabels([str(i + 1) for i in range(width)])
    ax.axhline(y=2, color='gray', linestyle='--', alpha=0.3)

    for p in range(width):
        ic = ic_per_pos[p]
        freqs = {n: freq_matrix[n][p] for n in NUCLEOTIDES}
        sorted_nucs = sorted(freqs, key=lambda n: freqs[n])

        y_offset = 0.0
        for nuc in sorted_nucs:
            letter_h = ic * freqs[nuc]
            if letter_h < 1e-4:
                continue
            # Draw a coloured rectangle
            rect = mpatches.FancyBboxPatch(
                (p + 0.05, y_offset), 0.9, letter_h,
                boxstyle="round,pad=0.02",
                facecolor=COLORS[nuc], edgecolor='white', linewidth=0.5
            )
            ax.add_patch(rect)
            # Text label
            if letter_h > 0.15:
                ax.text(p + 0.5, y_offset + letter_h / 2, nuc,
                        ha='center', va='center',
                        fontsize=max(6, min(20, int(letter_h * 18))),
                        fontweight='bold', color='white',
                        fontfamily='monospace')
            y_offset += letter_h

    # Legend
    legend_patches = [mpatches.Patch(color=COLORS[n], label=n)
                      for n in NUCLEOTIDES]
    ax.legend(handles=legend_patches, loc='upper right',
              ncol=4, fontsize=9, framealpha=0.8)

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {output_path}")
    plt.close()
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# 2. PWM Heatmap
# ──────────────────────────────────────────────────────────────────────────────

def plot_pwm_heatmap(pwm, title="PWM Heatmap", output_path=None,
                     figsize=(10, 3)):
    freq_matrix = pwm.frequency_matrix()
    data = np.array([[freq_matrix[n][p] for p in range(pwm.width)]
                     for n in NUCLEOTIDES])

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(data, cmap='Blues', aspect='auto', vmin=0, vmax=1)
    ax.set_yticks(range(4))
    ax.set_yticklabels(NUCLEOTIDES, fontsize=12, fontfamily='monospace')
    ax.set_xticks(range(pwm.width))
    ax.set_xticklabels([str(i + 1) for i in range(pwm.width)], fontsize=10)
    ax.set_xlabel("Position", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')

    for i in range(4):
        for j in range(pwm.width):
            val = data[i, j]
            color = 'white' if val > 0.6 else 'black'
            ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                    fontsize=8, color=color)

    plt.colorbar(im, ax=ax, label='Frequency')
    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {output_path}")
    plt.close()
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# 3. Gibbs Convergence
# ──────────────────────────────────────────────────────────────────────────────

def plot_gibbs_convergence(score_histories, labels=None, title="Gibbs Sampling Convergence",
                            output_path=None, figsize=(9, 4)):
    fig, ax = plt.subplots(figsize=figsize)
    colors_list = plt.cm.tab10(np.linspace(0, 1, len(score_histories)))

    for i, hist in enumerate(score_histories):
        label = labels[i] if labels else f"Run {i+1}"
        ax.plot(hist, color=colors_list[i], alpha=0.8, linewidth=1.5,
                label=label)

    ax.set_xlabel("Iteration", fontsize=12)
    ax.set_ylabel("Alignment Score", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(fontsize=9, ncol=min(5, len(score_histories)))
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {output_path}")
    plt.close()


# ──────────────────────────────────────────────────────────────────────────────
# 4. Algorithm Comparison
# ──────────────────────────────────────────────────────────────────────────────

def plot_algorithm_comparison(results, metric='information_content',
                               title=None, output_path=None, figsize=(8, 5)):
    """
    results: list of dicts with keys 'algorithm', 'dataset', and metric key.
    """
    if not results:
        return

    algorithms = list({r['algorithm'] for r in results})
    datasets   = list({r['dataset']   for r in results})

    x = np.arange(len(datasets))
    bar_w = 0.8 / max(1, len(algorithms))
    colors_list = plt.cm.Set2(np.linspace(0, 1, len(algorithms)))

    fig, ax = plt.subplots(figsize=figsize)

    for i, alg in enumerate(algorithms):
        vals = []
        for ds in datasets:
            match = [r[metric] for r in results
                     if r['algorithm'] == alg and r['dataset'] == ds]
            vals.append(match[0] if match else 0)

        offset = (i - len(algorithms) / 2 + 0.5) * bar_w
        bars = ax.bar(x + offset, vals, bar_w * 0.9, label=alg,
                      color=colors_list[i], edgecolor='white', linewidth=0.5)

        for bar, v in zip(bars, vals):
            if v > 0:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.02, f'{v:.2f}',
                        ha='center', va='bottom', fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(datasets, rotation=20, ha='right', fontsize=10)
    ax.set_ylabel(metric.replace('_', ' ').title(), fontsize=12)
    ax.set_title(title or f"Algorithm Comparison – {metric}", fontsize=14,
                 fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {output_path}")
    plt.close()


# ──────────────────────────────────────────────────────────────────────────────
# 5. Site Position Distribution
# ──────────────────────────────────────────────────────────────────────────────

def plot_site_distribution(sites, seq_length, title="Motif Site Distribution",
                            output_path=None, figsize=(9, 4)):
    positions = [s['position'] for s in sites]
    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(positions, bins=min(30, seq_length // 5),
            color='#3498db', edgecolor='white', alpha=0.85)
    ax.set_xlabel("Position in Sequence", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {output_path}")
    plt.close()


# ──────────────────────────────────────────────────────────────────────────────
# 6. Multi-motif IC comparison
# ──────────────────────────────────────────────────────────────────────────────

def plot_multi_motif_ic(motifs, title="Motif Information Content Comparison",
                         output_path=None, figsize=(8, 4)):
    """Plot IC per position for multiple motifs side by side."""
    n = len(motifs)
    fig, axes = plt.subplots(1, n, figsize=figsize, sharey=True)
    if n == 1:
        axes = [axes]

    for ax, m in zip(axes, motifs):
        pwm = m['pwm']
        freq_matrix = pwm.frequency_matrix()
        ic_per_pos = []
        for p in range(pwm.width):
            freqs = np.array([freq_matrix[nuc][p] for nuc in NUCLEOTIDES])
            freqs = np.clip(freqs, 1e-10, 1.0)
            ic = max(0, 2.0 + np.sum(freqs * np.log2(freqs)))
            ic_per_pos.append(ic)

        ax.bar(range(pwm.width), ic_per_pos, color='#8e44ad', alpha=0.8,
               edgecolor='white')
        ax.set_title(f"Motif {m['motif_id']}\n{m['consensus']}", fontsize=10)
        ax.set_xlabel("Position", fontsize=9)
        ax.set_ylim(0, 2.1)
        ax.axhline(y=2, color='gray', linestyle='--', alpha=0.3)

    axes[0].set_ylabel("IC (bits)", fontsize=11)
    fig.suptitle(title, fontsize=13, fontweight='bold')
    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {output_path}")
    plt.close()
