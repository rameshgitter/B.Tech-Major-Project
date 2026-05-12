#!/usr/bin/env python3
"""
pwm.py
Position Weight Matrix (PWM) / Position Frequency Matrix (PFM) implementation.
Includes scoring, consensus extraction, and information content calculation.

Author: Ramesh Chandra Soren (2022CSB086)
Supervisor: Dr. Surajeet Ghosh, IIEST Shibpur
"""

import math
import numpy as np
from collections import defaultdict

NUCLEOTIDES = ['A', 'C', 'G', 'T']
NUC_INDEX = {n: i for i, n in enumerate(NUCLEOTIDES)}


class PWM:
    """
    Position Weight Matrix for representing sequence motifs.
    Internally stores log-odds scores relative to background frequencies.
    """

    def __init__(self, width, pseudocount=0.1,
                 background=None):
        self.width = width
        self.pseudocount = pseudocount
        self.background = background or {'A': 0.25, 'C': 0.25,
                                         'G': 0.25, 'T': 0.25}
        # PFM: counts[position][nucleotide]
        self.counts = np.zeros((width, 4))
        self.nsites = 0
        # Will be computed after building
        self.log_odds = None
        self.info_content = None

    # ------------------------------------------------------------------ #
    #  Building the matrix
    # ------------------------------------------------------------------ #

    def add_site(self, site):
        """Add a single aligned motif occurrence (string of length=width)."""
        assert len(site) == self.width, \
            f"Site length {len(site)} != PWM width {self.width}"
        for pos, nuc in enumerate(site):
            if nuc in NUC_INDEX:
                self.counts[pos, NUC_INDEX[nuc]] += 1
        self.nsites += 1

    def add_sites(self, sites):
        for site in sites:
            self.add_site(site)

    def build_from_weighted(self, sequences, weights):
        """
        Build from fractional (EM-style) counts.
        sequences: list of strings of length width
        weights  : list of floats (probabilities)
        """
        self.counts = np.zeros((self.width, 4))
        self.nsites = 0
        for seq, w in zip(sequences, weights):
            for pos, nuc in enumerate(seq):
                if nuc in NUC_INDEX:
                    self.counts[pos, NUC_INDEX[nuc]] += w
            self.nsites += w
        self._compute_log_odds()

    def _compute_log_odds(self):
        """Convert counts to log-odds scores (with pseudocount smoothing)."""
        freq = (self.counts + self.pseudocount) / \
               (self.nsites + 4 * self.pseudocount)
        bg = np.array([self.background[n] for n in NUCLEOTIDES])
        with np.errstate(divide='ignore', invalid='ignore'):
            self.log_odds = np.log2(freq / bg)
        # Information content per position
        self.info_content = np.sum(freq * self.log_odds, axis=1)

    def finalize(self):
        self._compute_log_odds()

    # ------------------------------------------------------------------ #
    #  Scoring
    # ------------------------------------------------------------------ #

    def score(self, kmer):
        """Return log-odds score of a k-mer against this PWM."""
        if self.log_odds is None:
            self._compute_log_odds()
        s = 0.0
        for pos, nuc in enumerate(kmer):
            idx = NUC_INDEX.get(nuc)
            if idx is not None:
                s += self.log_odds[pos, idx]
            else:
                s += 0.0  # ambiguous base: no contribution
        return s

    def score_sequence(self, seq):
        """Score all windows in a sequence; return list of (position, score)."""
        results = []
        for i in range(len(seq) - self.width + 1):
            kmer = seq[i:i + self.width]
            results.append((i, self.score(kmer)))
        return results

    def best_match(self, seq):
        """Return (position, score, kmer) of the best-scoring window."""
        scored = self.score_sequence(seq)
        if not scored:
            return None
        pos, sc = max(scored, key=lambda x: x[1])
        return pos, sc, seq[pos:pos + self.width]

    # ------------------------------------------------------------------ #
    #  Summary
    # ------------------------------------------------------------------ #

    def consensus(self):
        """Return the consensus sequence (most probable nucleotide per position)."""
        if self.counts.sum() == 0:
            return 'N' * self.width
        return ''.join(NUCLEOTIDES[np.argmax(self.counts[p])]
                       for p in range(self.width))

    def total_information_content(self):
        if self.info_content is None:
            self._compute_log_odds()
        return float(np.sum(self.info_content))

    def frequency_matrix(self):
        """Return normalised frequency matrix as dict of lists."""
        total = self.counts + self.pseudocount
        total /= total.sum(axis=1, keepdims=True)
        return {n: list(total[:, NUC_INDEX[n]]) for n in NUCLEOTIDES}

    def __repr__(self):
        cons = self.consensus()
        ic = self.total_information_content() if self.info_content is not None else 0
        return f"PWM(width={self.width}, consensus='{cons}', IC={ic:.2f} bits)"

    # ------------------------------------------------------------------ #
    #  Serialisation
    # ------------------------------------------------------------------ #

    def to_dict(self):
        fm = self.frequency_matrix()
        return {
            'width': self.width,
            'nsites': int(self.nsites),
            'consensus': self.consensus(),
            'total_ic': self.total_information_content(),
            'frequency_matrix': fm
        }

    @classmethod
    def from_sites(cls, sites, pseudocount=0.1, background=None):
        """Convenience constructor: build a PWM from a list of aligned k-mers."""
        width = len(sites[0])
        pwm = cls(width, pseudocount=pseudocount, background=background)
        pwm.add_sites(sites)
        pwm.finalize()
        return pwm
