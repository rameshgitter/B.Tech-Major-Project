import math

NUCLEOTIDES = ['A','C','G','T']
IDX = {n:i for i,n in enumerate(NUCLEOTIDES)}


def build_pwm_from_motifs(motifs, alpha=1.0):
    """Build PWM probabilities (4 x k) with pseudocount alpha.
    motifs: list of equal-length strings
    returns: probs (list of 4 lists length k) and counts (list of 4 lists)
    """
    k = len(motifs[0])
    counts = [[0.0 for _ in range(k)] for _ in range(4)]
    for m in motifs:
        for j, b in enumerate(m):
            counts[IDX[b]][j] += 1.0

    probs = [[0.0 for _ in range(k)] for _ in range(4)]
    for j in range(k):
        col_sum = sum(counts[i][j] for i in range(4))
        denom = col_sum + 4.0 * alpha
        for i in range(4):
            probs[i][j] = (counts[i][j] + alpha) / denom
    return probs, counts


def pwm_log_prob(kmer, probs):
    """Return log probability of a kmer under the PWM (natural log)."""
    lp = 0.0
    for j, b in enumerate(kmer):
        lp += math.log(probs[IDX[b]][j])
    return lp


def consensus_from_pwm(probs):
    # probs is 4 x k
    k = len(probs[0])
    consensus = []
    for j in range(k):
        best_i = 0
        best_p = probs[0][j]
        for i in range(1, 4):
            if probs[i][j] > best_p:
                best_p = probs[i][j]
                best_i = i
        consensus.append(NUCLEOTIDES[best_i])
    return ''.join(consensus)
