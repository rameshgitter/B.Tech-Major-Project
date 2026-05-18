motif_demo README

This small demo runs a greedy motif search on `motifs.fasta` and writes results.

Files created by the demo:

- motif_demo_results.json: JSON with best k and motif details
- motif_demo_results.txt: human-readable summary
- motif_demo_pwm.png: PWM heatmap (requires matplotlib)
- motif_demo_logo.png: sequence logo (requires logomaker + matplotlib + pandas)

How to run:

1. (Optional) Create a virtual environment:
   python3 -m venv .venv
   source .venv/bin/activate

2. Install dependencies:
   pip install -r requirements.txt

3. Run the demo:
   python3 motif_demo.py

If you cannot install packages, the script will still produce `motif_demo_results.json` and `.txt` but will skip plotting.

Notes:

- The script searches k from 5 to 8 (you can change in the script).
- Ensure `motifs.fasta` contains only DNA lines (A/C/G/T). The demo ignores non-DNA lines.
