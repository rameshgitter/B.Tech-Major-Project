import math
import random

NUCLEOTIDES = ['A','C','G','T']

def information_content(probs, bg_prob=0.25):
    # probs: list of 4 lists (4 x k)
    k = len(probs[0])
    ic = 0.0
    ic_per_pos = []
    for j in range(k):
        s = 0.0
        for i in range(4):
            p = probs[i][j]
            if p > 0:
                s += p * math.log2(p / bg_prob)
        ic_per_pos.append(s)
        ic += s
    return ic, ic_per_pos

def estimate_evalue(sequences, k, best_score, trials=100, iterations=100, restarts=3):
    """Estimate an empirical E-value by shuffling sequences and counting how often a score >= best_score occurs."""
    from .gibbs_sampler import gibbs_sampler
    count = 0
    for t in range(trials):
        shuffled = []
        for s in sequences:
            lst = list(s)
            random.shuffle(lst)
            shuffled.append(''.join(lst))
        res = gibbs_sampler(shuffled, k, iterations=iterations, restarts=restarts)
        if res['score'] >= best_score:
            count += 1
    return count / max(1, trials)
