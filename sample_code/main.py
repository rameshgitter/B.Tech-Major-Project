#!/usr/bin/env python3
from pathlib import Path
from src.preprocessing import read_fasta
from src.gibbs_sampler import gibbs_sampler
from src.utils import save_json, ensure_dir
from src.visualize import save_pwm_heatmap_from_probs, save_sequence_logo_from_counts
from src.analysis import information_content, estimate_evalue

import argparse

parser = argparse.ArgumentParser(description='Gibbs sampler motif discovery demo')
parser.add_argument('--data', type=str, default='data/sequences/promoters.fasta', help='input fasta')
parser.add_argument('--k', type=int, default=6)
parser.add_argument('--iterations', type=int, default=200)
parser.add_argument('--restarts', type=int, default=10)
parser.add_argument('--out', type=str, default='results')
args = parser.parse_args()

p = Path(args.data)
if not p.exists():
    # fallback to repository motifs.fasta at parent
    fallback = Path(__file__).resolve().parents[1] / 'motifs.fasta'
    if fallback.exists():
        p = fallback
    else:
        raise FileNotFoundError(f'Input fasta not found: {args.data}')

seqs = read_fasta(str(p))
print(f'Loaded {len(seqs)} sequences from {p} (using k={args.k})')
# filter sequences shorter than k
seqs = [s for s in seqs if len(s) >= args.k]
if not seqs:
    raise SystemExit('No sequences long enough for chosen k')

# If only one sequence was found, try a simple fallback parser: accept any line that
# looks like a DNA string (only A/C/G/T) as a separate sequence. This handles files
# like the provided motifs.fasta which may contain many DNA lines without FASTA headers.
if len(seqs) == 1:
    try:
        with open(p, 'r') as fh:
            alt = []
            import re
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if re.fullmatch(r'[ACGTacgt]+', line):
                    alt.append(line.upper())
        alt = [s for s in alt if len(s) >= args.k]
        if len(alt) >= 2:
            seqs = alt
            print(f'Parsed {len(seqs)} DNA-only lines from {p} as sequences')
    except Exception:
        pass

if len(seqs) < 2:
    raise SystemExit('Need at least 2 sequences to run Gibbs sampler. Provide multiple sequences in FASTA or one-per-line DNA file.')

best = gibbs_sampler(seqs, args.k, iterations=args.iterations, restarts=args.restarts)

outdir = Path(args.out)
outdir.mkdir(parents=True, exist_ok=True)

save_json(best, str(outdir / 'gibbs_result.json'))
with open(outdir / 'consensus.txt', 'w') as f:
    f.write(best['consensus'] + '\n')
with open(outdir / 'motifs.txt', 'w') as f:
    for m in best['motifs']:
        f.write(m + '\n')

print('Saved results to', outdir)

# Visualization
try:
    pwm_png = outdir / 'pwm.png'
    logo_png = outdir / 'logo.png'
    probs = best['pwm_probs']
    counts = best['pwm_counts']
    ok1 = save_pwm_heatmap_from_probs(probs, str(pwm_png))
    ok2 = save_sequence_logo_from_counts(counts, str(logo_png))
    if ok1:
        print('Saved PWM heatmap to', pwm_png)
    if ok2:
        print('Saved sequence logo to', logo_png)
except Exception as e:
    print('Visualization failed:', e)

# Analysis
try:
    ic, ic_per_pos = information_content(best['pwm_probs'])
    with open(outdir / 'information_content.txt', 'w') as f:
        f.write(f'Total IC: {ic:.4f}\n')
        f.write('Per-position IC:\n')
        f.write('\n'.join(f'{x:.4f}' for x in ic_per_pos))
    print('Saved information content')
    # E-value estimation (fast, small trials by default)
    e = estimate_evalue(seqs, args.k, best['score'], trials=50, iterations=50, restarts=2)
    with open(outdir / 'evalue.txt', 'w') as f:
        f.write(f'Empirical E-value (fraction of shuffled trials >= observed): {e:.4f}\n')
    print('Saved E-value estimate')
except Exception as e:
    print('Analysis failed:', e)
