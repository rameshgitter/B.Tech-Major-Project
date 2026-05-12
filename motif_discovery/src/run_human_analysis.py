#!/usr/bin/env python3
"""
run_human_analysis.py — Motif discovery on human vertebrate TF datasets.
Real JASPAR data: CTCF (ChIP-seq), E2F1 (HT-SELEX + ChIP-seq), HOXA10 (HT-SELEX)
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(__file__))

from preprocess     import read_fasta, filter_sequences, compute_background_frequencies
from meme_algorithm import MEME
from gibbs_sampling import GibbsMotifFinder
from visualize      import (plot_sequence_logo, plot_pwm_heatmap,
                             plot_gibbs_convergence, plot_algorithm_comparison,
                             plot_multi_motif_ic)

RESULTS = "../results/real"
os.makedirs(RESULTS, exist_ok=True)

DATASETS = [
    {"name":"Human CTCF (ChIP-seq, MA0139.2)",
     "label":"ctcf","fasta":"../data/raw/jaspar_MA0139.2_CTCF.fasta",
     "width":10,"expected":"CCAGGGGGCG",
     "desc":"Real ChIP-seq: 19M genome-wide binding sites"},
    {"name":"Human E2F1 (HT-SELEX, MA0024.3)",
     "label":"e2f1_selex","fasta":"../data/raw/jaspar_MA0024.3_E2F1.fasta",
     "width":8,"expected":"TTTGGCGC",
     "desc":"Real HT-SELEX: cell-cycle regulatory factor"},
    {"name":"Human E2F1 (ChIP-seq, MA0024.2)",
     "label":"e2f1_chip","fasta":"../data/raw/jaspar_MA0024.2_E2F1_chipseq.fasta",
     "width":8,"expected":"GCGGGCGG",
     "desc":"Real ENCODE ChIP-seq"},
    {"name":"Human HOXA10 (HT-SELEX, MA0899.1)",
     "label":"hoxa10","fasta":"../data/raw/jaspar_MA0899.1_HOXA10.fasta",
     "width":8,"expected":"GTAATAAA",
     "desc":"Real HT-SELEX: HOX homeodomain TAAT core"},
    {"name":"All Human Vertebrate TFs (combined)",
     "label":"human_all","fasta":"../data/raw/real_human_vertebrate_TFs.fasta",
     "width":8,"expected":"GGCGGGCG",
     "desc":"CTCF + E2F1 (×2) + HOXA10 — 120 real sequences"},
]

print("\n" + "="*62)
print("  Human Vertebrate TF Motif Discovery")
print("  Real data: JASPAR CORE (jaspar.elixir.no) — April 2026")
print("="*62)

results = []
for cfg in DATASETS:
    raw = read_fasta(cfg["fasta"])
    filtered, _ = filter_sequences(raw, min_length=cfg["width"]+4)
    seqs = [s for _, s in filtered]
    bg   = compute_background_frequencies(filtered)

    print(f"\n  {cfg['name']}")
    print(f"  {cfg['desc']}  ({len(seqs)} seqs, width={cfg['width']})")
    print(f"  Expected motif: {cfg['expected']}")

    row = {"dataset":cfg["name"],"label":cfg["label"],
           "n_seqs":len(seqs),"width":cfg["width"],
           "expected":cfg["expected"]}

    # MEME
    t0 = time.time()
    meme = MEME(seqs, width=cfg["width"], n_motifs=2, model='ZOOPS',
                n_seeds=5, max_iter=15, background=bg)
    motifs = meme.fit(verbose=False)
    t = round(time.time()-t0,2)
    if motifs:
        m = motifs[0]
        print(f"  MEME  → {m['consensus']:<16} IC={m['information_content']:.2f}  "
              f"E={m['evalue']:.1e}  {t}s")
        row.update({"meme_consensus":m["consensus"],
                    "meme_ic":round(m["information_content"],3),
                    "meme_evalue":m["evalue"],"meme_time":t})
        plot_sequence_logo(m["pwm"],
            title=f"MEME — {cfg['name']}  |  {m['consensus']}",
            output_path=f"{RESULTS}/logo_meme_{cfg['label']}.png")
        plot_pwm_heatmap(m["pwm"],
            title=f"PWM Heatmap — {cfg['name']}",
            output_path=f"{RESULTS}/heatmap_{cfg['label']}.png")

    # Gibbs
    t0 = time.time()
    gibbs = GibbsMotifFinder(seqs, width=cfg["width"], n_runs=3,
                              n_iter=60, background=bg)
    gpwm, gscore = gibbs.fit(verbose=False)
    t = round(time.time()-t0,2)
    if gpwm:
        gc = gpwm.consensus(); gi = round(gpwm.total_information_content(),3)
        print(f"  Gibbs → {gc:<16} IC={gi:.2f}  "
              f"score={gscore:.2f}  {t}s")
        row.update({"gibbs_consensus":gc,"gibbs_ic":gi,
                    "gibbs_score":round(gscore,4),"gibbs_time":t})
        plot_sequence_logo(gpwm,
            title=f"Gibbs — {cfg['name']}  |  {gc}",
            output_path=f"{RESULTS}/logo_gibbs_{cfg['label']}.png")
        if gibbs.run_histories:
            plot_gibbs_convergence(gibbs.run_histories,
                labels=[f"Run {i+1}" for i in range(len(gibbs.run_histories))],
                title=f"Gibbs Convergence — {cfg['name']}",
                output_path=f"{RESULTS}/convergence_{cfg['label']}.png")

    results.append(row)

# Overall IC comparison for all real datasets
print("\n" + "="*62)
print("  Summary — Human Vertebrate TFs")
print(f"  {'Dataset':<38} {'MEME IC':>8}  {'Gibbs IC':>8}")
print("  " + "─"*56)
for r in results:
    mc = r.get("meme_ic",0); gc = r.get("gibbs_ic",0)
    print(f"  {r['dataset'][:38]:<38} {mc:>8.2f}  {gc:>8.2f}")

# Save combined IC comparison
all_results_path = f"{RESULTS}/results_real.json"
try:
    with open(all_results_path) as f:
        existing = json.load(f)
except:
    existing = []
existing.extend(results)
with open(all_results_path,"w") as f:
    json.dump(existing, f, indent=2)

# Master comparison plot
ic_cmp = []
for r in existing:
    ic_cmp.append({"algorithm":"MEME","dataset":r.get("label",r["dataset"][:12]),
                   "information_content":r.get("meme_ic",0)})
    ic_cmp.append({"algorithm":"Gibbs","dataset":r.get("label",r["dataset"][:12]),
                   "information_content":r.get("gibbs_ic",0)})

plot_algorithm_comparison(ic_cmp, metric="information_content",
    title="IC: MEME vs Gibbs — All Real JASPAR Datasets",
    output_path=f"{RESULTS}/comparison_ALL_ic.png", figsize=(14,5))

print(f"\n  All plots → {RESULTS}/")
print(f"  JSON     → {all_results_path}")
print("  Done.\n")
