import math

NUCLEOTIDES = ['A','C','G','T']
IDX = {n:i for i,n in enumerate(NUCLEOTIDES)}


def score_kmer_log(kmer, probs, bg_prob=0.25):
    """Compute log-likelihood ratio score = sum_j log( P(bj|PWM) / bg_prob ).
    Uses natural log.
    probs: list of 4 lists (4 x k)
    """
    s = 0.0
    for j, b in enumerate(kmer):
        p = probs[IDX[b]][j]
        s += math.log(p / bg_prob)
    return s


def total_log_likelihood(sequences, positions, k, probs, bg_prob=0.25):
    total = 0.0
    for seq, pos in zip(sequences, positions):
        kmer = seq[pos:pos+k]
        total += score_kmer_log(kmer, probs, bg_prob=bg_prob)
    return total
