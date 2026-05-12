#!/usr/bin/env python3
"""
run_all.py — Single entry point: fetches real data, runs all analyses,
             generates all plots and the final HTML report.

Usage:
    python run_all.py              # full pipeline
    python run_all.py --skip-fetch # skip NCBI/JASPAR fetch (use cached data)

Author: Ramesh Chandra Soren (2022CSB086) | IIEST Shibpur
"""
import os, sys, argparse, subprocess, time

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(ROOT, "data", "raw")
SRC  = os.path.join(ROOT, "src")

def run(script, cwd=None):
    cwd = cwd or ROOT
    print(f"\n{'─'*60}")
    print(f"  Running: {os.path.basename(script)}")
    print(f"{'─'*60}")
    t0 = time.time()
    result = subprocess.run([sys.executable, script], cwd=cwd, check=False)
    elapsed = time.time() - t0
    status = "OK" if result.returncode == 0 else "FAILED"
    print(f"  [{status}] {os.path.basename(script)} — {elapsed:.1f}s")
    return result.returncode == 0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-fetch", action="store_true",
                        help="Skip internet fetch, use cached data")
    args = parser.parse_args()

    print("\n" + "="*60)
    print("  BIOLOGICAL SEQUENCE MOTIF DISCOVERY")
    print("  B.Tech Major Project — IIEST Shibpur 2026")
    print("  Ramesh Chandra Soren (2022CSB086)")
    print("="*60)

    os.makedirs(os.path.join(ROOT, "results", "real"), exist_ok=True)
    os.makedirs(os.path.join(ROOT, "data", "processed"), exist_ok=True)

    # Step 1: Generate synthetic benchmarks
    run(os.path.join(DATA, "generate_synthetic_data.py"), cwd=DATA)

    # Step 2: Fetch / build real data
    if not args.skip_fetch:
        print("\n  [STEP 2] Building JASPAR real datasets...")
        run(os.path.join(DATA, "build_all_real_datasets.py"), cwd=DATA)
        run(os.path.join(DATA, "expand_real_data.py"), cwd=DATA)
        print("\n  [STEP 2b] Fetching NCBI + extra JASPAR (requires internet)...")
        run(os.path.join(DATA, "fetch_real_data.py"), cwd=DATA)
    else:
        print("\n  [STEP 2] Skipping fetch — using cached data")

    # Step 3: Run analyses
    run(os.path.join(SRC, "run_real_analysis.py"), cwd=SRC)
    run(os.path.join(SRC, "run_human_analysis.py"), cwd=SRC)
    run(os.path.join(SRC, "run_analysis.py"), cwd=SRC)

    # Step 4: Generate HTML report
    run(os.path.join(ROOT, "generate_report.py"), cwd=ROOT)

    print("\n" + "="*60)
    print("  PIPELINE COMPLETE")
    print(f"  Results  → {os.path.join(ROOT, 'results')}/")
    print(f"  Report   → {os.path.join(ROOT, 'report', 'index.html')}")
    print("="*60 + "\n")

if __name__ == "__main__":
    main()
