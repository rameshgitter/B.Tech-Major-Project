#!/usr/bin/env python3
"""
run_real_analysis.py
Full motif discovery pipeline on REAL biological datasets from JASPAR.
Runs MEME + Gibbs Sampling, evaluates results, generates all plots.

Author: Ramesh Chandra Soren (2022CSB086) — IIEST Shibpur
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(__file__))

from preprocess   import read_fasta, filter_sequences, compute_background_frequencies
from meme_algorithm import MEME
from gibbs_sampling import GibbsMotifFinder
from visualize    import (plot_sequence_logo, plot_pwm_heatmap,
                          plot_gibbs_convergence, plot_algorithm_comparison,
                          plot_site_distribution, plot_multi_motif_ic)

RESULTS = "../results/real"
os.makedirs(RESULTS, exist_ok=True)

# ── Real JASPAR datasets ──────────────────────────────────────────────────────
DATASETS = [
    {
        "name":   "Yeast TFs (S. cerevisiae)",
        "label":  "yeast",
        "fasta":  "../data/raw/real_yeast_TFs.fasta",
        "width":  8,
        "desc":   "ABF1/ABF2 PBM/DIP-chip (JASPAR MA0265.1, MA0266.1/2)",
        "expected_motif": "TCGTATAA",  # ABF1 core
    },
    {
        "name":   "Plant bZIP TFs (Arabidopsis)",
        "label":  "plant_bzip",
        "fasta":  "../data/raw/real_plant_bZIP.fasta",
        "width":  8,
        "desc":   "ABF1-4, ABI5 SELEX/PBM/ChIP/DAP-seq (JASPAR)",
        "expected_motif": "CACGTG",    # G-box / ACGT core
    },
    {
        "name":   "Plant B3/MADS TFs (Arabidopsis)",
        "label":  "plant_b3mads",
        "fasta":  "../data/raw/real_plant_B3_MADS.fasta",
        "width":  6,
        "desc":   "AGL3, ABI3, ABR1 SELEX/PBM/DAP-seq (JASPAR)",
        "expected_motif": "GCATGC",    # ABI3 core
    },
    {
        "name":   "All JASPAR TFs (combined)",
        "label":  "all_jaspar",
        "fasta":  "../data/raw/real_all_jaspar.fasta",
        "width":  8,
        "desc":   "16 JASPAR CORE matrices, all organisms combined",
        "expected_motif": "CACGTG",
    },
    # Also run on synthetic benchmarks for comparison
    {
        "name":   "Synthetic TATA-box (benchmark)",
        "label":  "synth_tata",
        "fasta":  "../data/raw/synthetic_tata_benchmark.fasta",
        "width":  8,
        "desc":   "Planted TATATAAG, 10% mutation rate",
        "expected_motif": "TATATAAG",
    },
    {
        "name":   "Synthetic E-box (benchmark)",
        "label":  "synth_ebox",
        "fasta":  "../data/raw/synthetic_ebox_benchmark.fasta",
        "width":  6,
        "desc":   "Planted CACGTG, 5% mutation rate",
        "expected_motif": "CACGTG",
    },
]

def run_dataset(cfg, verbose=True):
    label = cfg["label"]
    width = cfg["width"]
    name  = cfg["name"]

    # Load + preprocess
    raw = read_fasta(cfg["fasta"])
    filtered, _ = filter_sequences(raw, min_length=width + 4)
    seqs = [s for _, s in filtered]
    bg   = compute_background_frequencies(filtered)

    if not seqs:
        print(f"  [SKIP] No sequences after filtering for {name}")
        return None

    print(f"\n{'─'*62}")
    print(f"  {name}")
    print(f"  {cfg['desc']}")
    print(f"  {len(seqs)} sequences, width={width}")
    print(f"{'─'*62}")

    row = {"dataset": name, "label": label, "n_seqs": len(seqs), "width": width}

    # ── MEME ─────────────────────────────────────────────────────────────────
    t0 = time.time()
    meme = MEME(seqs, width=width, n_motifs=2, model='ZOOPS',
                n_seeds=5, max_iter=15, background=bg)
    meme_motifs = meme.fit(verbose=False)
    meme_time   = round(time.time() - t0, 2)

    if meme_motifs:
        m0 = meme_motifs[0]
        row.update({"meme_consensus": m0["consensus"],
                    "meme_ic": round(m0["information_content"], 3),
                    "meme_evalue": m0["evalue"],
                    "meme_time": meme_time})
        print(f"  MEME   → {m0['consensus']:20s}  IC={m0['information_content']:.2f}  "
              f"E={m0['evalue']:.1e}  t={meme_time}s")
    else:
        row.update({"meme_consensus": "–", "meme_ic": 0, "meme_evalue": 1, "meme_time": meme_time})
        print("  MEME   → no motif found")

    # ── Gibbs ─────────────────────────────────────────────────────────────────
    t0 = time.time()
    gibbs = GibbsMotifFinder(seqs, width=width, n_runs=3,
                              n_iter=60, background=bg)
    gibbs_pwm, gibbs_score = gibbs.fit(verbose=False)
    gibbs_time = round(time.time() - t0, 2)

    if gibbs_pwm:
        g_cons = gibbs_pwm.consensus()
        g_ic   = round(gibbs_pwm.total_information_content(), 3)
        row.update({"gibbs_consensus": g_cons, "gibbs_ic": g_ic,
                    "gibbs_score": round(gibbs_score, 4), "gibbs_time": gibbs_time})
        print(f"  Gibbs  → {g_cons:20s}  IC={g_ic:.2f}  "
              f"score={gibbs_score:.3f}  t={gibbs_time}s")
    else:
        row.update({"gibbs_consensus": "–", "gibbs_ic": 0,
                    "gibbs_score": 0, "gibbs_time": gibbs_time})
        print("  Gibbs  → no motif found")

    # ── Plots ─────────────────────────────────────────────────────────────────
    if meme_motifs:
        plot_sequence_logo(
            meme_motifs[0]["pwm"],
            title=f"MEME — {name}\nconsensus: {meme_motifs[0]['consensus']}",
            output_path=f"{RESULTS}/logo_meme_{label}.png")
        plot_pwm_heatmap(
            meme_motifs[0]["pwm"],
            title=f"PWM Heatmap — MEME — {name}",
            output_path=f"{RESULTS}/heatmap_{label}.png")
        if len(meme_motifs) > 1:
            plot_multi_motif_ic(
                meme_motifs,
                title=f"Motif IC Comparison — {name}",
                output_path=f"{RESULTS}/multi_ic_{label}.png")
        meme_sites = meme_motifs[0]["sites"]
        for i, s in enumerate(meme_sites):
            s.setdefault("seq_index", i)
        if meme_sites:
            avg_len = sum(len(s) for s in seqs) / len(seqs)
            plot_site_distribution(
                meme_sites, int(avg_len),
                title=f"Site Distribution — MEME — {name}",
                output_path=f"{RESULTS}/sites_{label}.png")

    if gibbs_pwm:
        plot_sequence_logo(
            gibbs_pwm,
            title=f"Gibbs Sampling — {name}\nconsensus: {gibbs_pwm.consensus()}",
            output_path=f"{RESULTS}/logo_gibbs_{label}.png")
        if gibbs.run_histories:
            plot_gibbs_convergence(
                gibbs.run_histories,
                labels=[f"Run {i+1}" for i in range(len(gibbs.run_histories))],
                title=f"Gibbs Convergence — {name}",
                output_path=f"{RESULTS}/convergence_{label}.png")

    return row


# ── Main ─────────────────────────────────────────────────────────────────────
print("\n" + "=" * 62)
print("  Motif Discovery on Real JASPAR Datasets")
print("  B.Tech Major Project | Ramesh Chandra Soren | IIEST Shibpur")
print("  Real data: JASPAR REST API (jaspar.elixir.no)")
print("=" * 62)

all_rows = []
for cfg in DATASETS:
    row = run_dataset(cfg)
    if row:
        all_rows.append(row)

# ── Summary table ─────────────────────────────────────────────────────────────
print("\n" + "=" * 62)
print(f"  {'Dataset':<32} {'MEME':^20} {'Gibbs':^20}")
print(f"  {'':32} {'Consensus':^12} {'IC':^6}  {'Consensus':^12} {'IC':^6}")
print("  " + "─" * 58)
for r in all_rows:
    mc = r.get("meme_consensus",  "–")[:12]
    mi = r.get("meme_ic",          0)
    gc = r.get("gibbs_consensus", "–")[:12]
    gi = r.get("gibbs_ic",         0)
    print(f"  {r['dataset'][:32]:<32} {mc:^12} {mi:^6.2f}  {gc:^12} {gi:^6.2f}")
print("=" * 62)

# Comparison plot (IC)
ic_cmp = []
for r in all_rows:
    ic_cmp.append({"algorithm":"MEME",  "dataset":r["label"],
                   "information_content": r.get("meme_ic", 0)})
    ic_cmp.append({"algorithm":"Gibbs", "dataset":r["label"],
                   "information_content": r.get("gibbs_ic", 0)})

plot_algorithm_comparison(
    ic_cmp, metric="information_content",
    title="Information Content: MEME vs Gibbs — Real JASPAR Datasets",
    output_path=f"{RESULTS}/comparison_real_ic.png")

# Save JSON
with open(f"{RESULTS}/results_real.json", "w") as f:
    json.dump(all_rows, f, indent=2)

print(f"\n  All plots → {RESULTS}/")
print(f"  Results  → {RESULTS}/results_real.json")
print("\n  Done.\n")
