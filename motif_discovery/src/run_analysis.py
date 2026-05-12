#!/usr/bin/env python3
"""
run_analysis.py
Master script: runs MEME and Gibbs Sampling on all datasets,
evaluates performance, and saves all results + plots.

Usage:
    cd src
    python run_analysis.py

Author: Ramesh Chandra Soren (2022CSB086)
Supervisor: Dr. Surajeet Ghosh, IIEST Shibpur
"""

import os
import sys
import json
import time

sys.path.insert(0, os.path.dirname(__file__))

from preprocess       import read_fasta, filter_sequences, compute_background_frequencies
from meme_algorithm   import MEME
from gibbs_sampling   import GibbsMotifFinder
from evaluate         import (site_metrics, read_true_positions_from_fasta,
                               motif_similarity_score, print_evaluation_table)
from visualize        import (plot_sequence_logo, plot_pwm_heatmap,
                               plot_gibbs_convergence,
                               plot_algorithm_comparison,
                               plot_site_distribution,
                               plot_multi_motif_ic)

# ──────────────────────────────────────────────────────────────────────────────
# Config
# ──────────────────────────────────────────────────────────────────────────────

DATASETS = [
    {
        'name':          'TATA-box (Yeast)',
        'fasta':         '../data/raw/synthetic_tata_benchmark.fasta',
        'true_width':    8,
        'true_motif':    'TATATAAG',
        'has_ground_truth': True,
    },
    {
        'name':          'E-box (MYC)',
        'fasta':         '../data/raw/synthetic_ebox_benchmark.fasta',
        'true_width':    6,
        'true_motif':    'CACGTG',
        'has_ground_truth': True,
    },
    {
        'name':          'NF-kB',
        'fasta':         '../data/raw/synthetic_nfkb_benchmark.fasta',
        'true_width':    10,
        'true_motif':    'GGGACTTTCC',
        'has_ground_truth': True,
    },
    {
        'name':          'Zinc Finger',
        'fasta':         '../data/raw/synthetic_znf_benchmark.fasta',
        'true_width':    8,
        'true_motif':    'GCGTGGCG',
        'has_ground_truth': True,
    },
]

RESULTS_DIR = '../results'
os.makedirs(RESULTS_DIR, exist_ok=True)

# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

def load_dataset(cfg):
    raw = read_fasta(cfg['fasta'])
    filtered, _ = filter_sequences(raw, min_length=cfg['true_width'] + 10)
    sequences = [s for _, s in filtered]
    bg = compute_background_frequencies(filtered)
    return sequences, bg

# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main():
    all_eval_results = []
    comparison_data  = []
    meme_gibbs_comparison = []

    print("\n" + "=" * 70)
    print("  Biological Sequence Motif Discovery – Full Analysis Pipeline")
    print("  B.Tech Major Project | IIEST Shibpur | 2026")
    print("=" * 70)

    for cfg in DATASETS:
        name    = cfg['name']
        width   = cfg['true_width']
        print(f"\n{'─'*70}")
        print(f"  Dataset: {name}  |  Width: {width}  |  True: {cfg['true_motif']}")
        print(f"{'─'*70}")

        sequences, bg = load_dataset(cfg)
        print(f"  Loaded {len(sequences)} sequences")

        true_pos = {}
        if cfg['has_ground_truth']:
            true_pos = read_true_positions_from_fasta(cfg['fasta'])

        # ────────────────────────────────────────────────────────────────────
        # MEME
        # ────────────────────────────────────────────────────────────────────
        print(f"\n  [MEME]")
        t0 = time.time()
        meme = MEME(sequences, width=width, n_motifs=1, model='ZOOPS',
                    n_seeds=4, max_iter=15, background=bg)
        meme_motifs = meme.fit(verbose=False)
        meme_time = time.time() - t0

        meme_best = meme_motifs[0] if meme_motifs else None
        meme_sites = meme_best['sites'] if meme_best else []
        meme_ic    = meme_best['information_content'] if meme_best else 0

        # Attach seq_index to MEME sites
        for i, site in enumerate(meme_sites):
            if 'seq_index' not in site:
                site['seq_index'] = i

        print(f"  Consensus : {meme_best['consensus'] if meme_best else 'N/A'}")
        print(f"  IC        : {meme_ic:.2f} bits")
        print(f"  Time      : {meme_time:.1f}s")

        # ────────────────────────────────────────────────────────────────────
        # Gibbs Sampling
        # ────────────────────────────────────────────────────────────────────
        print(f"\n  [Gibbs Sampling]")
        t0 = time.time()
        gibbs = GibbsMotifFinder(sequences, width=width, n_runs=3,
                                  n_iter=50, background=bg)
        gibbs_pwm, gibbs_score = gibbs.fit(verbose=False)
        gibbs_time = time.time() - t0
        gibbs_sites = gibbs.best_sites
        gibbs_ic = gibbs_pwm.total_information_content() if gibbs_pwm else 0

        print(f"  Consensus : {gibbs_pwm.consensus() if gibbs_pwm else 'N/A'}")
        print(f"  IC        : {gibbs_ic:.2f} bits")
        print(f"  Time      : {gibbs_time:.1f}s")

        # ────────────────────────────────────────────────────────────────────
        # Evaluation
        # ────────────────────────────────────────────────────────────────────
        if true_pos:
            meme_eval  = site_metrics(meme_sites, true_pos, width)
            gibbs_eval = site_metrics(gibbs_sites, true_pos, width)

            print(f"\n  Evaluation:")
            print(f"  {'Metric':<20} {'MEME':>10} {'Gibbs':>10}")
            print(f"  {'─'*40}")
            for k in ['sensitivity', 'ppv', 'f1', 'nCC', 'site_accuracy']:
                print(f"  {k:<20} {meme_eval[k]:>10.4f} {gibbs_eval[k]:>10.4f}")

            all_eval_results.append({
                'algorithm': 'MEME', 'dataset': name,
                'information_content': round(meme_ic, 3),
                **meme_eval
            })
            all_eval_results.append({
                'algorithm': 'Gibbs', 'dataset': name,
                'information_content': round(gibbs_ic, 3),
                **gibbs_eval
            })

        comparison_data.append({
            'dataset':   name,
            'MEME':      {'ic': meme_ic,  'time': meme_time,
                          'consensus': meme_best['consensus'] if meme_best else ''},
            'Gibbs':     {'ic': gibbs_ic, 'time': gibbs_time,
                          'consensus': gibbs_pwm.consensus() if gibbs_pwm else ''},
        })

        # ────────────────────────────────────────────────────────────────────
        # Save plots
        # ────────────────────────────────────────────────────────────────────
        safe_name = name.replace(' ', '_').replace('(', '').replace(')', '')

        if meme_best:
            plot_sequence_logo(
                meme_best['pwm'],
                title=f"MEME – {name}  |  consensus: {meme_best['consensus']}",
                output_path=f"{RESULTS_DIR}/logo_meme_{safe_name}.png"
            )
            plot_pwm_heatmap(
                meme_best['pwm'],
                title=f"PWM Heatmap – MEME – {name}",
                output_path=f"{RESULTS_DIR}/heatmap_meme_{safe_name}.png"
            )
            if len(meme_motifs) > 1:
                plot_multi_motif_ic(
                    meme_motifs,
                    title=f"MEME Motif IC – {name}",
                    output_path=f"{RESULTS_DIR}/multi_ic_meme_{safe_name}.png"
                )

        if gibbs_pwm:
            plot_sequence_logo(
                gibbs_pwm,
                title=f"Gibbs – {name}  |  consensus: {gibbs_pwm.consensus()}",
                output_path=f"{RESULTS_DIR}/logo_gibbs_{safe_name}.png"
            )
            # Use score history from already-run gibbs
            histories = gibbs.run_histories if gibbs.run_histories else [[0]]

            plot_gibbs_convergence(
                histories,
                labels=[f'Run {i+1}' for i in range(len(histories))],
                title=f"Gibbs Convergence – {name}",
                output_path=f"{RESULTS_DIR}/convergence_gibbs_{safe_name}.png"
            )

        if meme_sites:
            avg_seq_len = sum(len(s) for s in sequences) / len(sequences)
            plot_site_distribution(
                meme_sites, int(avg_seq_len),
                title=f"MEME Site Distribution – {name}",
                output_path=f"{RESULTS_DIR}/sites_meme_{safe_name}.png"
            )

    # ──────────────────────────────────────────────────────────────────────
    # Summary comparison plot
    # ──────────────────────────────────────────────────────────────────────
    print(f"\n{'='*70}")
    print("  Overall Evaluation Table")
    print_evaluation_table(all_eval_results)

    ic_results = []
    for row in comparison_data:
        ic_results.append({
            'algorithm': 'MEME',  'dataset': row['dataset'],
            'information_content': row['MEME']['ic']
        })
        ic_results.append({
            'algorithm': 'Gibbs', 'dataset': row['dataset'],
            'information_content': row['Gibbs']['ic']
        })

    plot_algorithm_comparison(
        ic_results, metric='information_content',
        title="Information Content: MEME vs Gibbs Sampling",
        output_path=f"{RESULTS_DIR}/comparison_ic.png"
    )

    # F1 comparison
    if all_eval_results:
        plot_algorithm_comparison(
            all_eval_results, metric='f1',
            title="F1 Score: MEME vs Gibbs Sampling",
            output_path=f"{RESULTS_DIR}/comparison_f1.png"
        )

    # Save JSON results
    json_path = f"{RESULTS_DIR}/evaluation_results.json"
    with open(json_path, 'w') as f:
        json.dump({'evaluation': all_eval_results,
                   'comparison': comparison_data}, f, indent=2)
    print(f"\n  Results saved to {json_path}")
    print(f"  All plots saved in {RESULTS_DIR}/")
    print("\n  Analysis complete.\n")


if __name__ == "__main__":
    main()
