import math
import json
import re
from pathlib import Path

NUCLEOTIDES = ['A','C','G','T']
IDX = {n:i for i,n in enumerate(NUCLEOTIDES)}


def parse_pfm(path):
    """Parse a simple JASPAR-like PFM.
    Accepts formats like:
      A [ 12 3 4 ... ]
      C [ ... ]
      G [ ... ]
      T [ ... ]
    or whitespace-separated rows starting with A/C/G/T.
    Returns probs (4 x k) as lists of floats.
    """
    rows = { 'A':[], 'C':[], 'G':[], 'T':[] }
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                continue
            m = re.match(r'^([ACGTacgt])[^0-9-]*(.*)$', line)
            if not m:
                continue
            base = m.group(1).upper()
            nums = re.findall(r'[-+]?[0-9]*\.?[0-9]+', m.group(2))
            if not nums:
                continue
            rows[base] = [float(x) for x in nums]
    # validate
    lengths = {len(v) for v in rows.values() if v}
    if not lengths:
        raise ValueError('No numeric rows found in PFM')
    if len(lengths) != 1:
        raise ValueError('PFM rows have inconsistent lengths: %r' % lengths)
    k = lengths.pop()
    counts = [rows['A'], rows['C'], rows['G'], rows['T']]
    # convert to probabilities per column
    probs = [[0.0]*k for _ in range(4)]
    for j in range(k):
        colsum = sum(counts[i][j] for i in range(4))
        if colsum <= 0:
            for i in range(4):
                probs[i][j] = 0.25
        else:
            for i in range(4):
                probs[i][j] = counts[i][j] / colsum
    return probs


def load_discovered_pwm(results_json_path):
    with open(results_json_path) as f:
        data = json.load(f)
    # expects 'pwm_probs' key with 4 x k lists
    probs = data.get('pwm_probs')
    if probs is None:
        raise ValueError('No pwm_probs found in results JSON')
    # Ensure order is [A,C,G,T]
    return probs


def flatten_pwm_segment(probs, start, length):
    flat = []
    for j in range(start, start+length):
        for i in range(4):
            flat.append(probs[i][j])
    return flat


def pearson(x, y):
    if len(x) != len(y):
        raise ValueError('length mismatch')
    n = len(x)
    mx = sum(x)/n
    my = sum(y)/n
    num = sum((x[i]-mx)*(y[i]-my) for i in range(n))
    denx = math.sqrt(sum((x[i]-mx)**2 for i in range(n)))
    deny = math.sqrt(sum((y[i]-my)**2 for i in range(n)))
    if denx == 0 or deny == 0:
        return 0.0
    return num / (denx * deny)


def best_alignment_similarity(found_probs, pfm_probs):
    # found_probs: 4 x kf, pfm_probs: 4 x kp
    kf = len(found_probs[0])
    kp = len(pfm_probs[0])
    best = { 'score': -2.0, 'offset': None, 'overlap': 0 }
    # slide pfm across found (allow pfm to be shorter or longer)
    for offset in range(-kp+1, kf):
        # overlap region in found: start_f .. end_f-1
        start_f = max(0, offset)
        start_p = max(0, -offset)
        overlap = min(kf - start_f, kp - start_p)
        if overlap <= 0:
            continue
        xf = flatten_pwm_segment(found_probs, start_f, overlap)
        xp = flatten_pwm_segment(pfm_probs, start_p, overlap)
        score = pearson(xf, xp)
        if score > best['score']:
            best = { 'score': score, 'offset': offset, 'overlap': overlap }
    return best


def main(results_json, pfm_path, out_path=None):
    found = load_discovered_pwm(results_json)
    pfm = parse_pfm(pfm_path)
    res = best_alignment_similarity(found, pfm)
    out = {
        'results_json': str(results_json),
        'pfm_path': str(pfm_path),
        'best_score': res['score'],
        'offset': res['offset'],
        'overlap': res['overlap']
    }
    if out_path:
        with open(out_path, 'w') as f:
            import json as _json
            _json.dump(out, f, indent=2)
    else:
        print('Best similarity:', res['score'])
        print('Offset:', res['offset'])
        print('Overlap columns:', res['overlap'])
    return out


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--results', default='motif_project/results_part2/gibbs_result.json')
    p.add_argument('--pfm')
    p.add_argument('--out')
    args = p.parse_args()
    if not args.pfm:
        raise SystemExit('Please provide --pfm path to JASPAR PFM file')
    main(args.results, args.pfm, args.out)
