# Biological Sequence Motif Discovery
### B.Tech Major Project | IIEST Shibpur | 2026

**Student:** Ramesh Chandra Soren (2022CSB086)
**Supervisor:** Dr. Surajeet Ghosh
**Department:** Computer Science and Technology

---

## Quick Start

```bash
pip install numpy matplotlib requests biopython
python run_all.py                 # full pipeline
python run_all.py --skip-fetch    # skip internet fetch
open report/index.html            # view HTML report
```

---

## Project Structure

```
motif_discovery/
├── run_all.py                        Single entry point for full pipeline
├── generate_report.py                Self-contained HTML report generator
├── requirements.txt
│
├── data/raw/
│   ├── fetch_real_data.py            Fetch JASPAR API + NCBI Entrez
│   ├── build_all_real_datasets.py    Build FASTA from fetched PFMs
│   ├── expand_real_data.py           Add human TFs (CTCF, E2F1, HOXA10)
│   ├── generate_synthetic_data.py    Benchmark datasets with planted motifs
│   ├── real_human_vertebrate_TFs.fasta  CTCF + E2F1 + HOXA10 (JASPAR CORE)
│   ├── real_yeast_TFs.fasta          S. cerevisiae ABF1/ABF2
│   ├── real_plant_bZIP.fasta         Arabidopsis ABF1-4, ABI5
│   ├── real_plant_B3_MADS.fasta      Arabidopsis AGL3, ABI3, ABR1
│   ├── real_all_jaspar.fasta         All 22 JASPAR matrices combined
│   ├── jaspar_*.fasta                Per-TF FASTA files (22 files)
│   ├── synthetic_*.fasta             Synthetic benchmarks (4 files)
│   └── jaspar_real_metadata.json     Full provenance for all 22 matrices
│
├── src/
│   ├── pwm.py                        Position Weight Matrix
│   ├── preprocess.py                 FASTA I/O and cleaning
│   ├── meme_algorithm.py             Full MEME EM (OOPS/ZOOPS, multi-motif)
│   ├── gibbs_sampling.py             Gibbs Sampler (multi-start)
│   ├── evaluate.py                   Sensitivity, PPV, F1, nCC
│   ├── visualize.py                  Logos, heatmaps, convergence plots
│   ├── run_real_analysis.py          Real JASPAR data pipeline
│   ├── run_human_analysis.py         Human vertebrate TF pipeline
│   └── run_analysis.py               Synthetic benchmark pipeline
│
├── results/real/                     58 plots + results JSON
├── report/index.html                 Self-contained HTML report (3.7 MB)
├── docs/RESULTS_REPORT.md            Detailed findings
└── notebooks/motif_discovery_demo.ipynb
```

---

## Real Data — JASPAR CORE (22 matrices, live API fetch)

| Matrix ID  | TF      | Organism         | Experiment   | PubMed   |
|------------|---------|------------------|--------------|----------|
| MA0139.2   | CTCF    | Homo sapiens     | ChIP-seq     | 21106759 |
| MA0024.3   | E2F1    | Homo sapiens     | HT-SELEX     | 9372931  |
| MA0024.2   | E2F1    | Homo sapiens     | ChIP-seq     | 17908821 |
| MA0899.1   | HOXA10  | Homo sapiens     | HT-SELEX     | 18585359 |
| MA0265.1   | ABF1    | S. cerevisiae    | PBM/DIP-chip | 19111667 |
| MA0266.1/2 | ABF2    | S. cerevisiae    | PBM          | 19111667 |
| MA0001.1   | AGL3    | A. thaliana      | SELEX        | 7632923  |
| MA0570.1   | ABF1    | A. thaliana      | SELEX        | 10636868 |
| MA0930.1-3 | ABF3    | A. thaliana      | PBM+ChIP-seq | 11005831 |
| MA0931.1/2 | ABI5    | A. thaliana      | PBM          | 12376636 |
| MA0941.1   | ABF2    | A. thaliana      | PBM          | 11005831 |
| MA1659.1/2 | ABF4    | A. thaliana      | DAP-seq      | 25220462 |
| MA0564.1/2 | ABI3    | A. thaliana      | PBM          | 30183137 |
| MA1244.1   | ABR1    | A. thaliana      | DAP-seq      | 9756931  |
| MA0020.2   | Dof2    | Zea mays         | SELEX        | 10074718 |
| MA0021.1   | Dof3    | Zea mays         | SELEX        | 10074718 |

---

## Results Summary

| Dataset               | MEME consensus | IC    | Gibbs consensus | IC     | E-value    |
|-----------------------|----------------|-------|-----------------|--------|------------|
| CTCF (ChIP-seq)       | CAGGGGGCAC     | 9.00  | CACCAGGGGG      | 11.29  | 4.7e-138   |
| E2F1 (HT-SELEX)       | TTGGCGCC       | 10.11 | TTGGCGCC        | 11.12  | 5.7e-152   |
| E2F1 (ChIP-seq)       | GGCGGGAG       | 8.47  | GGCGGGAG        | 9.47   | 1.7e-130   |
| HOXA10 (HT-SELEX)     | GTAATAAA       | 6.52  | GTAATAAA        | 7.38   | 6.7e-101   |
| Plant bZIP G-box      | ACACGTGT       | 8.64  | ACACGTGT        | 8.85   | 1.0e-300   |
| Plant B3/MADS         | GCATGC         | 4.69  | GCATGC          | 5.15   | 3.8e-190   |
| Yeast ABF1/ABF2       | CTCTAGAA       | 6.63  | CTAGATAT        | 6.72   | 1.2e-191   |

Gibbs Sampling outperformed MEME in IC across all 11 datasets (avg +1.2 bits).

---

## References

1. Castro-Mondragon et al. (2022). JASPAR 2022. Nucleic Acids Res 50:D165.
2. Baydar et al. (2026). JASPAR 2026. Nucleic Acids Res 54:D184.
3. Bailey & Elkan (1994). MEME. ISMB 2:28-36.
4. Lawrence et al. (1993). Gibbs Sampling. Science 262:208-214.
5. Tompa et al. (2005). Benchmarking motif discovery. Nat Biotechnol 23:137.
