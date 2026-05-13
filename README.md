# 🧬 Biological Sequence Motif Discovery

A bioinformatics project implementing multiple algorithms for discovering recurring patterns (motifs) in DNA, RNA, and protein sequences.

**Department of Computer Science and Technology, IIEST Shibpur**
Under the supervision of **Dr. Surajeet Ghosh**

---

## 📌 What are Biological Motifs?

Biological motifs are short, recurring patterns (typically 6–20 nucleotides or amino acids) in biological sequences that are functionally significant — such as transcription factor binding sites, promoter elements, and protein active sites.

---

## 🔬 Algorithms Implemented

| Algorithm | Type | Description |
|---|---|---|
| **Gibbs Sampling** | Probabilistic | Iteratively samples motif positions using probability distributions; generates Position Weight Matrices (PWM) |
| **Expectation Maximization (MEME)** | Probabilistic | Randomly initializes a PWM, then alternates E-step and M-step until convergence |
| **Markov Models** | Statistical | Models nucleotide occurrence probability based on preceding positions |
| **Needleman-Wunsch** | Exact / Dynamic Programming | Global sequence alignment with configurable match/mismatch/gap scoring |
| **Suffix Trees** | Exact | Represents all suffixes for fast pattern matching; finds motif positions in O(n) |
| **Box-Links** | Exact | Divides sequences into segments with inter-segment pointers for faster search |
| **Random Projections** | Heuristic | Maps sequences to lower-dimensional space using Gaussian projection matrices |
| **Uniform Projections** | Heuristic | Improves on random projections using a discrepancy threshold to reduce search space |

---

## 🚀 Getting Started

```bash
git clone https://github.com/<your-username>/bio-motif-discovery.git
cd bio-motif-discovery
pip install -r requirements.txt
```

---

## 📁 Project Structure

```
bio-motif-discovery/
├── data/               # Sample DNA/RNA sequences (FASTA format)
├── algorithms/         # Individual algorithm implementations
├── utils/              # PWM generation, sequence encoding, visualization
├── results/            # Output motifs and sequence logos
└── main.py             # Entry point
```

---

## 🧪 Example

Input sequence: `AGTCGTCATGA`

- **Markov Model** predicts next nucleotide based on transition probabilities
- **Suffix Tree** finds pattern `GT` at positions `[2, 5]`
- **Needleman-Wunsch** aligns `AGTCA` vs `GTCA` → score: `3`

---

## 📊 Motif Representation

- **Consensus Sequence** — most frequent nucleotide per position (e.g., `TATAAT`)
- **Position Weight Matrix (PWM)** — frequency table per position
- **Sequence Logo** — visual information-content plot

---

## 🌐 Applications

- Gene regulation & transcription factor binding site discovery
- Drug target identification
- Cancer genomics
- Personalized medicine & comparative genomics

---

## 📚 References

- He et al. (2021). *A survey on deep learning in DNA/RNA motif mining.* Briefings in Bioinformatics.
- Marschall & Rahmann (2009). *Efficient exact motif discovery.* Bioinformatics.
- Federico & Pisanti (2009). *Suffix tree characterization of maximal motifs.* Theoretical Computer Science.
- Pradhan (2008). *Motif discovery in biological sequences.*
