#!/usr/bin/env python3
"""
gibbs_sampling.py
Full implementation of the Gibbs Sampling algorithm for motif discovery.

Reference:
  Lawrence et al. (1993). Detecting subtle sequence signals: a Gibbs sampling
  strategy for multiple alignment. Science, 262(5131):208-214.

Author: Ramesh Chandra Soren (2022CSB086)
Supervisor: Dr. Surajeet Ghosh, IIEST Shibpur
"""

import math
import random
import numpy as np
from copy import deepcopy
from pwm import PWM, NUCLEOTIDES, NUC_INDEX

random.seed(42)
np.random.seed(42)


class GibbsSampler:
    """
    Gibbs Sampling for motif discovery.

    Algorithm:
      1. Randomly initialise one start position per sequence.
      2. For each sequence in turn:
         a. Remove it from the current alignment.
         b. Re-estimate the motif model from the remaining sequences.
         c. Score every window in the held-out sequence.
         d. Sample a new start position proportionally to the scores.
      3. Track the best alignment seen across all iterations.

    Parameters
    ----------
    sequences     : list of strings
    width         : motif width (k)
    background    : dict {A/C/G/T: frequency}
    pseudocount   : Laplace smoothing for PWM
    n_iterations  : Gibbs sweep iterations
    """

    def __init__(self, sequences, width=8, background=None,
                 pseudocount=0.1, n_iterations=200):
        self.sequences = sequences
        self.width = width
        self.background = background or {n: 0.25 for n in NUCLEOTIDES}
        self.pseudocount = pseudocount
        self.n_iterations = n_iterations

        # Validate: sequences must be at least width long
        self.valid_idx = [i for i, s in enumerate(sequences)
                          if len(s) >= width]
        self.n_valid = len(self.valid_idx)

        # Current alignment: start positions for each valid sequence
        self.positions = {}
        self.best_positions = {}
        self.best_score = -np.inf
        self.best_pwm = None

        self.score_history = []

    # ------------------------------------------------------------------ #
    #  Initialisation
    # ------------------------------------------------------------------ #

    def _random_init(self):
        for idx in self.valid_idx:
            seq = self.sequences[idx]
            max_start = len(seq) - self.width
            self.positions[idx] = random.randint(0, max_start)

    # ------------------------------------------------------------------ #
    #  Build PWM from all sequences except one
    # ------------------------------------------------------------------ #

    def _build_pwm_excluding(self, exclude_idx):
        counts = np.full((self.width, 4), self.pseudocount)
        nsites = 0

        for idx in self.valid_idx:
            if idx == exclude_idx:
                continue
            seq = self.sequences[idx]
            pos = self.positions[idx]
            kmer = seq[pos:pos + self.width]
            for p, nuc in enumerate(kmer):
                j = NUC_INDEX.get(nuc)
                if j is not None:
                    counts[p, j] += 1
            nsites += 1

        # Normalise
        freq = counts / (counts.sum(axis=1, keepdims=True))
        return freq   # shape (width, 4)

    # ------------------------------------------------------------------ #
    #  Score a k-mer against a frequency matrix
    # ------------------------------------------------------------------ #

    def _kmer_log_odds(self, kmer, freq):
        score = 0.0
        for p, nuc in enumerate(kmer):
            j = NUC_INDEX.get(nuc)
            if j is not None:
                motif_p = max(freq[p, j], 1e-300)
                bg_p = max(self.background.get(nuc, 0.25), 1e-300)
                score += math.log(motif_p / bg_p)
        return score

    # ------------------------------------------------------------------ #
    #  Sample a new position for sequence `idx`
    # ------------------------------------------------------------------ #

    def _sample_position(self, idx, freq):
        seq = self.sequences[idx]
        n_pos = len(seq) - self.width + 1
        if n_pos <= 0:
            return 0

        log_scores = np.array([
            self._kmer_log_odds(seq[i:i + self.width], freq)
            for i in range(n_pos)
        ])

        # Numerically stable softmax
        log_scores -= log_scores.max()
        probs = np.exp(log_scores)
        probs /= probs.sum()

        return np.random.choice(n_pos, p=probs)

    # ------------------------------------------------------------------ #
    #  Alignment score (sum of best-window log-odds)
    # ------------------------------------------------------------------ #

    def _alignment_score(self):
        freq = self._build_pwm_for_all()
        total = 0.0
        for idx in self.valid_idx:
            seq = self.sequences[idx]
            pos = self.positions[idx]
            kmer = seq[pos:pos + self.width]
            total += self._kmer_log_odds(kmer, freq)
        return total

    def _build_pwm_for_all(self):
        counts = np.full((self.width, 4), self.pseudocount)
        for idx in self.valid_idx:
            seq = self.sequences[idx]
            pos = self.positions[idx]
            kmer = seq[pos:pos + self.width]
            for p, nuc in enumerate(kmer):
                j = NUC_INDEX.get(nuc)
                if j is not None:
                    counts[p, j] += 1
        return counts / counts.sum(axis=1, keepdims=True)

    # ------------------------------------------------------------------ #
    #  Main sampling loop
    # ------------------------------------------------------------------ #

    def run(self):
        self._random_init()

        for iteration in range(self.n_iterations):
            # Randomly shuffle the order in which sequences are updated
            order = list(self.valid_idx)
            random.shuffle(order)

            for idx in order:
                freq = self._build_pwm_excluding(idx)
                new_pos = self._sample_position(idx, freq)
                self.positions[idx] = new_pos

            # Evaluate current alignment
            score = self._alignment_score()
            self.score_history.append(score)

            if score > self.best_score:
                self.best_score = score
                self.best_positions = deepcopy(self.positions)

        # Build final PWM from best alignment
        self.positions = self.best_positions
        freq = self._build_pwm_for_all()
        # Convert to PWM object
        self.best_pwm = PWM(self.width, self.pseudocount, self.background)
        self.best_pwm.counts = (
            freq * (self.n_valid + 4 * self.pseudocount)
        )
        self.best_pwm.nsites = self.n_valid
        self.best_pwm.finalize()

        return self.best_pwm, self.best_score

    # ------------------------------------------------------------------ #
    #  Retrieve results
    # ------------------------------------------------------------------ #

    def get_sites(self):
        """Return list of dicts {seq_idx, position, kmer, score}."""
        sites = []
        if not self.best_positions:
            return sites
        freq = self._build_pwm_for_all()
        for idx in self.valid_idx:
            seq = self.sequences[idx]
            pos = self.best_positions[idx]
            kmer = seq[pos:pos + self.width]
            score = self._kmer_log_odds(kmer, freq)
            sites.append({
                'seq_index': idx,
                'position': pos,
                'kmer': kmer,
                'score': score
            })
        return sorted(sites, key=lambda x: -x['score'])


# ──────────────────────────────────────────────────────────────────────────────
# Multi-start wrapper
# ──────────────────────────────────────────────────────────────────────────────

class GibbsMotifFinder:
    """
    Runs Gibbs Sampler multiple times (multi-start) to escape local optima.

    Parameters
    ----------
    sequences  : list of strings
    width      : motif width
    n_runs     : number of independent Gibbs restarts
    n_iter     : Gibbs iterations per run
    background : nucleotide background frequencies
    """

    def __init__(self, sequences, width=8, n_runs=5, n_iter=200,
                 background=None, pseudocount=0.1):
        self.sequences = sequences
        self.width = width
        self.n_runs = n_runs
        self.n_iter = n_iter
        self.background = background or {n: 0.25 for n in NUCLEOTIDES}
        self.pseudocount = pseudocount
        self.best_pwm = None
        self.best_score = -np.inf
        self.best_sites = []
        self.run_scores = []
        self.run_histories = []  # score_history per run

    def fit(self, verbose=True):
        if verbose:
            print(f"\n[Gibbs] width={self.width}, "
                  f"n_runs={self.n_runs}, n_iter={self.n_iter}")

        for run in range(self.n_runs):
            sampler = GibbsSampler(
                self.sequences, self.width, self.background,
                self.pseudocount, self.n_iter
            )
            pwm, score = sampler.run()
            self.run_scores.append(score)
            self.run_histories.append(sampler.score_history)

            if verbose:
                print(f"  Run {run + 1:2d}/{self.n_runs}  "
                      f"score={score:.3f}  "
                      f"consensus={pwm.consensus()}")

            if score > self.best_score:
                self.best_score = score
                self.best_pwm = deepcopy(pwm)
                self.best_sites = sampler.get_sites()
                self.best_sampler = sampler

        if verbose:
            print(f"\nBest consensus : {self.best_pwm.consensus()}")
            print(f"Best score     : {self.best_score:.3f}")
            print(f"IC             : "
                  f"{self.best_pwm.total_information_content():.2f} bits")
        return self.best_pwm, self.best_score

    def summary(self):
        if self.best_pwm is None:
            print("No motif found yet. Run .fit() first.")
            return
        print("\n" + "=" * 60)
        print("Gibbs Sampling Results")
        print("=" * 60)
        print(f"Consensus        : {self.best_pwm.consensus()}")
        print(f"Best alignment   : {self.best_score:.4f}")
        print(f"Information cont.: "
              f"{self.best_pwm.total_information_content():.2f} bits")
        print(f"Sites discovered : {len(self.best_sites)}")
        print("-" * 40)
        print("Top 5 sites:")
        for s in self.best_sites[:5]:
            print(f"  seq={s['seq_index']:3d}  "
                  f"pos={s['position']:4d}  "
                  f"kmer={s['kmer']}  "
                  f"score={s['score']:.3f}")
        print("=" * 60)


# ──────────────────────────────────────────────────────────────────────────────
# CLI entry point
# ──────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys
    sys.path.insert(0, '.')
    from preprocess import read_fasta, filter_sequences, compute_background_frequencies

    fasta = "../data/raw/synthetic_ebox_benchmark.fasta"
    raw = read_fasta(fasta)
    seqs_data, _ = filter_sequences(raw, min_length=20)
    sequences = [s for _, s in seqs_data]
    bg = compute_background_frequencies(seqs_data)

    print(f"Running Gibbs Sampling on {len(sequences)} sequences...")
    finder = GibbsMotifFinder(sequences, width=6, n_runs=5, n_iter=300,
                               background=bg)
    finder.fit(verbose=True)
    finder.summary()
