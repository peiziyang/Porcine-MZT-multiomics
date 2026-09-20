"""SRA FASTQ download + salmon quantification pipeline.
Downloads a sample's paired FASTQ.gz from EBI, runs salmon quant, cleans up.
Usage: python sra_salmon_pipeline.py [SRR_id] [condition] [stage]
"""
import os, sys, subprocess, time

SAMPLE = sys.argv[1] if len(sys.argv) > 1 else "SRR13435642"
CONDITION = sys.argv[2] if len(sys.argv) > 2 else "?"
STAGE = sys.argv[3] if len(sys.argv) > 3 else "?"

TMPDIR = "/mnt/e/Workbuddy/2026-07-27-11-58-27/data/raw/sra_tmp"
OUTDIR = "/mnt/e/Workbuddy/2026-07-27-11-58-27/data/processed/salmon_out"
INDEX  = "/mnt/e/Workbuddy/2026-07-27-11-58-27/data/raw/salmon_index"
WSL = "wsl --distribution Ubuntu --exec bash -c"

os.makedirs(TMPDIR.replace("/mnt/e/", "E:/"), exist_ok=True)
os.makedirs(OUTDIR.replace("/mnt/e/", "E:/"), exist_ok=True)

start = time.time()
print(f"[{SAMPLE}] {CONDITION} {STAGE} - start", flush=True)

# Step 1: Download R1 and R2 FASTQ from EBI using WSL+curl
# EBI URL pattern: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR134/042/SRR13435642/SRR13435642_{1,2}.fastq.gz
prefix = f"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR134/042/{SAMPLE}/{SAMPLE}"

for read_num in ["1", "2"]:
    dest = f"{TMPDIR}/{SAMPLE}_{read_num}.fastq.gz"
    if os.path.exists(dest):
        print(f"  [{SAMPLE}] R{read_num} already exists, skipping", flush=True)
        continue
    url = f"{prefix}_{read_num}.fastq.gz"
    print(f"  [{SAMPLE}] Downloading R{read_num}...", flush=True)
    # Use wsl's curl with explicit full path (no shell interpolation issues)
    cmd = f'{WSL} "curl -sS -L -o {dest} {url}"'
    rc = os.system(cmd)
    if rc != 0:
        print(f"  [{SAMPLE}] curl R{read_num} FAILED (rc={rc})", flush=True)
        # Try wget as fallback
        cmd2 = f'{WSL} "wget -q -O {dest} {url}"'
        rc2 = os.system(cmd2)
        if rc2 != 0:
            print(f"  [{SAMPLE}] wget R{read_num} also FAILED", flush=True)
            sys.exit(1)

# Step 2: Run salmon quant
quant_dir = f"{OUTDIR}/{SAMPLE}_{CONDITION}_{STAGE}"
r1 = f"{TMPDIR}/{SAMPLE}_1.fastq.gz"
r2 = f"{TMPDIR}/{SAMPLE}_2.fastq.gz"

if os.path.exists(f"{quant_dir}/quant.sf"):
    print(f"  [{SAMPLE}] Already quantified, skipping", flush=True)
else:
    print(f"  [{SAMPLE}] Running salmon quant...", flush=True)
    cmd = f'{WSL} "salmon quant -i {INDEX} -l A -1 {r1} -2 {r2} --validateMappings -p 4 -o {quant_dir}"'
    rc = os.system(cmd)
    if rc != 0:
        print(f"  [{SAMPLE}] salmon quant FAILED (rc={rc})", flush=True)
        sys.exit(1)

# Step 3: Clean up
for read_num in ["1", "2"]:
    p = f"{TMPDIR}/{SAMPLE}_{read_num}.fastq.gz"
    if os.path.exists(p):
        os.remove(p)

elapsed = time.time() - start
print(f"  [{SAMPLE}] DONE in {elapsed/60:.1f} min", flush=True)
