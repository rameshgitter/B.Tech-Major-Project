import random
import math
from .pwm import build_pwm_from_motifs, consensus_from_pwm
from .scoring import total_log_likelihood

random.seed(0)


def gibbs_sampler(sequences, k, iterations=200, restarts=10, alpha=1.0, bg_prob=0.25):
    """Run Gibbs sampling multiple restarts; return best result dict.
    sequences: list of strings (ACGT)
    returns: dict with keys: positions, motifs, pwm_probs, pwm_counts, consensus, score, k
    """
    best = None
    n = len(sequences)
    for r in range(restarts):
        # initialize random positions
        positions = [random.randint(0, len(seq)-k) for seq in sequences]
        for it in range(iterations):
            for i in range(n):
                # build motifs excluding i
                motifs = [sequences[j][positions[j]:positions[j]+k] for j in range(n) if j != i]
                probs, counts = build_pwm_from_motifs(motifs, alpha=alpha)
                seq = sequences[i]
                L = len(seq) - k + 1
                logps = [0.0] * L
                for start in range(L):
                    kmer = seq[start:start+k]
                    # log probability under PWM
                    lp = 0.0
                    for j, b in enumerate(kmer):
                        idx = {'A':0, 'C':1, 'G':2, 'T':3}[b]
                        lp += math.log(probs[idx][j])
                    logps[start] = lp
                # stable softmax
                maxlp = max(logps)
                weights = [math.exp(x - maxlp) for x in logps]
                total_w = sum(weights)
                if total_w <= 0:
                    weights = [1.0] * len(weights)
                    total_w = len(weights)
                probs_norm = [w / total_w for w in weights]
                # sample new position
                rnum = random.random()
                cum = 0.0
                chosen = 0
                for idx, p in enumerate(probs_norm):
                    cum += p
                    if rnum <= cum:
                        chosen = idx
                        break
                positions[i] = chosen
        # after iterations compute final pwm and score
        motifs_final = [sequences[j][positions[j]:positions[j]+k] for j in range(n)]
        pwm_probs, pwm_counts = build_pwm_from_motifs(motifs_final, alpha=alpha)
        score = total_log_likelihood(sequences, positions, k, pwm_probs, bg_prob=bg_prob)
        if best is None or score > best['score']:
            best = {
                'positions': positions.copy(),
                'motifs': motifs_final,
                'pwm_probs': pwm_probs,
                'pwm_counts': pwm_counts,
                'consensus': consensus_from_pwm(pwm_probs),
                'score': float(score),
                'k': k,
                'restarts': r,
            }
    # convert lists to JSON-serializable
    best['pwm_probs'] = best['pwm_probs']
    best['pwm_counts'] = best['pwm_counts']
    return best
