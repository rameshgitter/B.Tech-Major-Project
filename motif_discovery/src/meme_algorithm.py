#!/usr/bin/env python3
"""
meme_algorithm.py
Full implementation of the MEME (Multiple EM for Motif Elicitation) algorithm.
Supports OOPS (One Occurrence Per Sequence) and ZOOPS (Zero or One Occurrence) models.

Reference:
  Bailey & Elkan (1994). Fitting a mixture model by expectation maximization
  to discover motifs in biopolymers. ISMB, 2:28-36.

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


# ──────────────────────────────────────────────────────────────────────────────
# Helper: background model
# ──────────────────────────────────────────────────────────────────────────────

def _bg_score(kmer, background):
    s = 0.0
    for nuc in kmer:
        p = background.get(nuc, 0.25)
        s += math.log(p + 1e-300)
    return s


# ──────────────────────────────────────────────────────────────────────────────
# Single EM run
# ──────────────────────────────────────────────────────────────────────────────

class MEMERun:
    """
    One EM run starting from a single seed/candidate motif.
    """

    def __init__(self, sequences, width, background, model='ZOOPS',
                 pseudocount=0.1, max_iter=100, tol=1e-6):
        self.sequences = sequences          # list of strings
        self.width = width
        self.background = background
        self.model = model
        self.pseudocount = pseudocount
        self.max_iter = max_iter
        self.tol = tol
        self.pwm = None
        self.log_likelihood = -np.inf

    # ------------------------------------------------------------------ #
    #  Initialisation from a seed k-mer
    # ------------------------------------------------------------------ #

    def _init_from_seed(self, seed):
        """Bootstrap the PWM from a single k-mer."""
        self.pwm = PWM(self.width, self.pseudocount, self.background)
        self.pwm.add_site(seed)
        # Add small uniform counts to avoid zero probabilities
        self.pwm.counts += 0.1
        self.pwm.nsites += 1
        self.pwm.finalize()

    # ------------------------------------------------------------------ #
    #  E-step: compute posterior probabilities
    # ------------------------------------------------------------------ #

    def _e_step(self):
        """
        For each position i in each sequence s, compute
            z[s][i] = P(motif starts at i | sequence s, current PWM)
        Returns list-of-arrays and per-sequence total responsibility.
        """
        all_z = []
        for seq in self.sequences:
            L = len(seq)
            n_pos = L - self.width + 1
            if n_pos <= 0:
                all_z.append(np.zeros(1))
                continue

            raw = np.zeros(n_pos)
            for i in range(n_pos):
                kmer = seq[i:i + self.width]
                motif_ll = sum(
                    math.log(max(self.pwm.counts[p, NUC_INDEX.get(kmer[p], 0)]
                                 / (self.pwm.nsites + 4 * self.pseudocount),
                                 1e-300))
                    for p in range(self.width)
                    if kmer[p] in NUC_INDEX
                )
                bg_ll = _bg_score(kmer, self.background)
                raw[i] = motif_ll - bg_ll  # log-odds

            # Convert to probabilities (softmax)
            raw -= raw.max()
            probs = np.exp(raw)

            if self.model == 'OOPS':
                probs /= probs.sum() + 1e-300
            else:  # ZOOPS: motif may be absent
                # Lambda = prior probability a sequence contains the motif
                lam = 0.5
                probs = lam * probs / (probs.sum() + 1e-300)

            all_z.append(probs)
        return all_z

    # ------------------------------------------------------------------ #
    #  M-step: update PWM from responsibilities
    # ------------------------------------------------------------------ #

    def _m_step(self, all_z):
        new_counts = np.full((self.width, 4), self.pseudocount)
        total_weight = 0.0
        for seq, z in zip(self.sequences, all_z):
            n_pos = len(seq) - self.width + 1
            for i in range(min(n_pos, len(z))):
                kmer = seq[i:i + self.width]
                for p, nuc in enumerate(kmer):
                    idx = NUC_INDEX.get(nuc)
                    if idx is not None:
                        new_counts[p, idx] += z[i]
                total_weight += z[i]

        self.pwm.counts = new_counts
        self.pwm.nsites = total_weight + 4 * self.pseudocount
        self.pwm.finalize()

    # ------------------------------------------------------------------ #
    #  Log-likelihood of data under current model
    # ------------------------------------------------------------------ #

    def _compute_log_likelihood(self):
        ll = 0.0
        for seq in self.sequences:
            n_pos = len(seq) - self.width + 1
            if n_pos <= 0:
                continue
            site_lls = []
            for i in range(n_pos):
                kmer = seq[i:i + self.width]
                site_lls.append(self.pwm.score(kmer))
            ll += max(site_lls)   # ZOOPS approximation
        return ll

    # ------------------------------------------------------------------ #
    #  Main EM loop
    # ------------------------------------------------------------------ #

    def run(self, seed):
        self._init_from_seed(seed)
        prev_ll = -np.inf

        for iteration in range(self.max_iter):
            all_z = self._e_step()
            self._m_step(all_z)
            ll = self._compute_log_likelihood()

            if abs(ll - prev_ll) < self.tol:
                break
            prev_ll = ll

        self.log_likelihood = ll
        return self.pwm, ll, all_z


# ──────────────────────────────────────────────────────────────────────────────
# MEME main class (multiple seeds, multiple motifs)
# ──────────────────────────────────────────────────────────────────────────────

class MEME:
    """
    MEME algorithm: finds multiple motifs by iterating EM over candidate seeds.

    Parameters
    ----------
    sequences : list of str
    width     : int, motif width
    n_motifs  : int, number of motifs to find
    model     : 'OOPS' or 'ZOOPS'
    n_seeds   : number of random seed starts per motif
    max_iter  : maximum EM iterations per run
    """

    def __init__(self, sequences, width=8, n_motifs=3,
                 model='ZOOPS', n_seeds=10, max_iter=50,
                 background=None, pseudocount=0.1):
        self.sequences = sequences
        self.width = width
        self.n_motifs = n_motifs
        self.model = model
        self.n_seeds = n_seeds
        self.max_iter = max_iter
        self.background = background or self._estimate_background()
        self.pseudocount = pseudocount
        self.motifs = []            # list of (PWM, log_likelihood, sites)

    def _estimate_background(self):
        from collections import Counter
        c = Counter()
        for seq in self.sequences:
            c.update(seq)
        total = sum(c.values()) + 1e-10
        return {n: c.get(n, 1) / total for n in NUCLEOTIDES}

    # ------------------------------------------------------------------ #
    #  Seed generation
    # ------------------------------------------------------------------ #

    def _generate_seeds(self):
        """Extract candidate seed k-mers from all sequences."""
        seeds = []
        for seq in self.sequences:
            # Sample at most 2 seeds per sequence at random positions
            import random as _r
            n_pos = len(seq) - self.width + 1
            if n_pos <= 0:
                continue
            positions = _r.sample(range(n_pos), min(2, n_pos))
            for i in positions:
                seeds.append(seq[i:i + self.width])
        _r.shuffle(seeds)
        return seeds[:self.n_seeds]

    # ------------------------------------------------------------------ #
    #  Erasing: mask identified motif occurrences
    # ------------------------------------------------------------------ #

    def _erase_motif(self, pwm, sequences):
        """
        Replace best-scoring window in each sequence with 'N' characters
        to allow discovery of the next motif.
        """
        erased = []
        for seq in sequences:
            pos, sc, kmer = pwm.best_match(seq) or (0, 0, '')
            if pos is not None and sc > 0:
                masked = seq[:pos] + 'N' * self.width + seq[pos + self.width:]
                erased.append(masked)
            else:
                erased.append(seq)
        return erased

    # ------------------------------------------------------------------ #
    #  Collect best motif sites
    # ------------------------------------------------------------------ #

    def _collect_sites(self, pwm):
        sites = []
        for seq in self.sequences:
            result = pwm.best_match(seq)
            if result:
                pos, sc, kmer = result
                sites.append({'position': pos, 'score': sc, 'kmer': kmer})
        return sites

    # ------------------------------------------------------------------ #
    #  E-value estimation (simplified)
    # ------------------------------------------------------------------ #

    def _evalue(self, ll, width, n_sequences):
        """
        Simplified E-value: lower is better.
        Real MEME uses a more elaborate formula; this is a reasonable proxy.
        """
        # Binomial approximation
        expected = n_sequences * 50 / (4 ** width)  # very rough
        evalue = max(1e-300, expected * math.exp(-ll))
        return evalue

    # ------------------------------------------------------------------ #
    #  Main fit method
    # ------------------------------------------------------------------ #

    def fit(self, verbose=True):
        working_seqs = list(self.sequences)

        for motif_idx in range(self.n_motifs):
            if verbose:
                print(f"\n[MEME] Searching for motif {motif_idx + 1} "
                      f"(width={self.width}, model={self.model})")

            seeds = self._generate_seeds()
            best_pwm = None
            best_ll = -np.inf
            best_z = None

            for seed_idx, seed in enumerate(seeds):
                if 'N' in seed:
                    continue
                runner = MEMERun(
                    working_seqs, self.width, self.background,
                    model=self.model, pseudocount=self.pseudocount,
                    max_iter=self.max_iter
                )
                try:
                    pwm, ll, z = runner.run(seed)
                except Exception:
                    continue

                if ll > best_ll:
                    best_ll = ll
                    best_pwm = deepcopy(pwm)
                    best_z = z

            if best_pwm is None:
                if verbose:
                    print(f"  No motif found for round {motif_idx + 1}")
                break

            sites = self._collect_sites(best_pwm)
            evalue = self._evalue(best_ll, self.width, len(working_seqs))
            ic = best_pwm.total_information_content()

            self.motifs.append({
                'motif_id': motif_idx + 1,
                'pwm': best_pwm,
                'log_likelihood': best_ll,
                'evalue': evalue,
                'information_content': ic,
                'consensus': best_pwm.consensus(),
                'sites': sites,
                'n_sites': len(sites)
            })

            if verbose:
                print(f"  Consensus: {best_pwm.consensus()}")
                print(f"  Log-likelihood: {best_ll:.4f}")
                print(f"  E-value: {evalue:.2e}")
                print(f"  Information content: {ic:.2f} bits")
                print(f"  Sites found: {len(sites)}")

            # Erase this motif before searching for the next
            working_seqs = self._erase_motif(best_pwm, working_seqs)

        return self.motifs

    def summary(self):
        print("\n" + "=" * 60)
        print("MEME Results Summary")
        print("=" * 60)
        for m in self.motifs:
            print(f"\nMotif {m['motif_id']}: {m['consensus']}")
            print(f"  IC={m['information_content']:.2f} bits  "
                  f"E-value={m['evalue']:.2e}  "
                  f"Sites={m['n_sites']}")
        print("=" * 60)


# ──────────────────────────────────────────────────────────────────────────────
# CLI entry point
# ──────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys
    sys.path.insert(0, '.')
    from preprocess import read_fasta, filter_sequences, compute_background_frequencies

    fasta = "../data/raw/synthetic_tata_benchmark.fasta"
    raw = read_fasta(fasta)
    seqs_data, _ = filter_sequences(raw, min_length=20)
    sequences = [s for _, s in seqs_data]
    bg = compute_background_frequencies(seqs_data)

    print(f"Running MEME on {len(sequences)} sequences...")
    meme = MEME(sequences, width=8, n_motifs=2, model='ZOOPS',
                n_seeds=15, background=bg)
    motifs = meme.fit(verbose=True)
    meme.summary()
