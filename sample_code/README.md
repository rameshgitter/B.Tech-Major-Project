Gibbs Sampling motif discovery demo

Structure:

- data/ (put jaspar PFMs and sequences here)
- src/ (implementation)
- results/ (output)

Quick run (create a venv recommended):

python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python3 main.py --data ../motifs.fasta --k 5 --iterations 200 --restarts 5 --out results

Notes:

- The demo uses only numpy; plotting is optional and not included here.
- Place your FASTA at data/sequences/promoters.fasta to use real data.

## PFM Comparison

To compare the discovered PWM with a JASPAR PFM, place the PFM file at:

    motif_project/data/jaspar/TP53.pfm

Then run (from project root, activate the venv first):

    python3 -m src.compare_pfm --results results_part2/gibbs_result.json --pfm motif_project/data/jaspar/TP53.pfm --out results_part2/pwm_vs_jaspar.json

The script returns a Pearson similarity score and best alignment offset.
