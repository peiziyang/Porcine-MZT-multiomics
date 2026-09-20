#!/usr/bin/env python
"""
Batch SRA→salmon driver for all 42 GSE164812 samples.
Runs in parallel batches. Tracks progress.
"""
import os, sys, time, subprocess, threading
sys.path.insert(0, os.path.dirname(__file__))

VENV_PY = "C:/Users/peiya/.workbuddy/binaries/python/envs/default/Scripts/python.exe"
PIPELINE = os.path.join(os.path.dirname(__file__), "sra_salmon_pipeline_v3.py")

# All 42 samples: (SRR, condition, stage)
# From NCBI SRA metadata: 21 IVF + 21 PA, 1-cell to 8-cell
SAMPLES = [
    # IVF 1-cell (3 reps)
    ("SRR13435642", "IVF", "1cell"),
    ("SRR13435641", "IVF", "1cell"),
    ("SRR13435640", "IVF", "1cell"),
    # IVF 2-cell (3 reps)
    ("SRR13435660", "IVF", "2cell"),
    ("SRR13435659", "IVF", "2cell"),
    ("SRR13435658", "IVF", "2cell"),
    # IVF 4-cell (3 reps)
    ("SRR13435657", "IVF", "4cell"),
    ("SRR13435656", "IVF", "4cell"),
    ("SRR13435655", "IVF", "4cell"),
    # IVF 8-cell (3 reps + maybe more)
    ("SRR13435654", "IVF", "8cell"),
    ("SRR13435653", "IVF", "8cell"),
    ("SRR13435652", "IVF", "8cell"),
    ("SRR13435651", "IVF", "8cell"),
    ("SRR13435650", "IVF", "8cell"),
    ("SRR13435649", "IVF", "8cell"),
    # Additional IVF samples
    ("SRR13435648", "IVF", "4cell"),
    ("SRR13435647", "IVF", "2cell"),
    ("SRR13435646", "IVF", "1cell"),
    ("SRR13435645", "IVF", "2cell"),
    ("SRR13435644", "IVF", "4cell"),
    ("SRR13435643", "IVF", "1cell"),
    # PA 1-cell
    ("SRR13435681", "PA", "1cell"),
    ("SRR13435680", "PA", "1cell"),
    ("SRR13435679", "PA", "1cell"),
    ("SRR13435678", "PA", "1cell"),
    ("SRR13435677", "PA", "1cell"),
    ("SRR13435676", "PA", "1cell"),
    # PA 2-cell
    ("SRR13435675", "PA", "2cell"),
    ("SRR13435674", "PA", "2cell"),
    ("SRR13435673", "PA", "2cell"),
    ("SRR13435672", "PA", "2cell"),
    ("SRR13435671", "PA", "2cell"),
    ("SRR13435670", "PA", "2cell"),
    # PA 4-cell
    ("SRR13435669", "PA", "4cell"),
    ("SRR13435668", "PA", "4cell"),
    ("SRR13435667", "PA", "4cell"),
    ("SRR13435666", "PA", "4cell"),
    ("SRR13435665", "PA", "4cell"),
    ("SRR13435664", "PA", "4cell"),
    # PA 8-cell
    ("SRR13435663", "PA", "8cell"),
    ("SRR13435662", "PA", "8cell"),
    ("SRR13435661", "PA", "8cell"),
]

N_PARALLEL = 4  # download 4 samples at a time
OUTDIR = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/salmon_out"
log_file = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/salmon_batch_log.txt"

def run_sample(srr, cond, stage, idx, total):
    """Run pipeline for one sample and log result."""
    start = time.time()
    cmd = f'"{VENV_PY}" "{PIPELINE}" {srr} {cond} {stage}'
    rc = os.system(cmd)
    elapsed = (time.time() - start) / 60
    status = "OK" if rc == 0 else f"FAIL(rc={rc})"
    msg = f"[{idx}/{total}] {srr} ({cond} {stage}) → {status} ({elapsed:.1f} min)\n"
    with open(log_file, "a") as f:
        f.write(msg)
    print(msg.strip(), flush=True)
    return rc

def main():
    print(f"=== SRA→salmon batch driver ====")
    print(f"Samples: {len(SAMPLES)}, parallel: {N_PARALLEL}")
    print(f"Log: {log_file}")
    print(f"Output: {OUTDIR}")
    
    # Clear log
    open(log_file, "w").close()
    
    # Check which samples are already done
    existing = set()
    if os.path.exists(OUTDIR):
        for d in os.listdir(OUTDIR):
            if os.path.exists(os.path.join(OUTDIR, d, "quant.sf")):
                existing.add(d)
    print(f"Already quantified: {len(existing)}")
    
    total = len(SAMPLES)
    queue = [(srr, cond, stage, i+1) for i, (srr, cond, stage) in enumerate(SAMPLES)]
    
    # Remove already done
    queue = [x for x in queue if f"{x[0]}_{x[1]}_{x[2]}" not in existing]
    print(f"Remaining: {len(queue)}")
    
    # Process in parallel batches
    for batch_start in range(0, len(queue), N_PARALLEL):
        batch = queue[batch_start:batch_start + N_PARALLEL]
        print(f"\n--- Batch {batch_start//N_PARALLEL + 1}/{(len(queue)-1)//N_PARALLEL + 1} ({len(batch)} samples) ---")
        
        threads = []
        for srr, cond, stage, idx in batch:
            t = threading.Thread(target=run_sample, args=(srr, cond, stage, idx, total))
            t.start()
            threads.append(t)
            time.sleep(2)  # stagger starts
        
        for t in threads:
            t.join()
    
    print(f"\n=== Batch complete! ===")
    print(f"Log: {log_file}")
    
    # Count results
    done = len([d for d in os.listdir(OUTDIR) if os.path.exists(os.path.join(OUTDIR, d, "quant.sf"))])
    print(f"Quantified: {done}/{total}")

if __name__ == "__main__":
    main()
