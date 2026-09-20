"""SRA FASTQ download + salmon quantification pipeline v2.
Uses Python urllib for reliable download + WSL salmon for quant.
Usage: python sra_salmon_pipeline.py [SRR_id] [condition] [stage]
"""
import os, sys, urllib.request, time, shutil

SAMPLE = sys.argv[1] if len(sys.argv) > 1 else "SRR13435642"
COND = sys.argv[2] if len(sys.argv) > 2 else "?"
STAGE = sys.argv[3] if len(sys.argv) > 3 else "?"

# Windows paths
TMPDIR_WIN = "E:/Workbuddy/2026-07-27-11-58-27/data/raw/sra_tmp"
OUTDIR_WIN = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/salmon_out"
# WSL paths
TMPDIR_WSL = "/mnt/e/Workbuddy/2026-07-27-11-58-27/data/raw/sra_tmp"
INDEX_WSL = "/mnt/e/Workbuddy/2026-07-27-11-58-27/data/raw/salmon_index"
OUTDIR_WSL = "/mnt/e/Workbuddy/2026-07-27-11-58-27/data/processed/salmon_out"

os.makedirs(TMPDIR_WIN, exist_ok=True)
os.makedirs(OUTDIR_WIN, exist_ok=True)

start = time.time()
print(f"[{SAMPLE}] {COND} {STAGE} - start", flush=True)

# Step 1: Download R1 and R2 FASTQ via Python urllib (reliable, 2.7 MB/s)
SIZE_TOTAL = 0
for read_num in ["1", "2"]:
    dest = f"{TMPDIR_WIN}/{SAMPLE}_{read_num}.fastq.gz"
    if os.path.exists(dest) and os.path.getsize(dest) > 10 * 1024 * 1024:
        print(f"  [{SAMPLE}] R{read_num} exists ({os.path.getsize(dest)/1e6:.0f} MB), skipping", flush=True)
        SIZE_TOTAL += os.path.getsize(dest)
        continue
    
    url = f"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR134/042/{SAMPLE}/{SAMPLE}_{read_num}.fastq.gz"
    print(f"  [{SAMPLE}] Downloading R{read_num} from EBI...", flush=True)
    
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        resp = urllib.request.urlopen(req, timeout=120)
        total = int(resp.headers.get("content-length", 0))
        print(f"  [{SAMPLE}] R{read_num} size: {total/1e9:.2f} GB", flush=True)
        
        CHUNK = 8 * 1024 * 1024
        downloaded = 0
        with open(dest, "wb") as f:
            while True:
                chunk = resp.read(CHUNK)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)
                if downloaded % (50 * 1024 * 1024) < CHUNK:
                    pct = downloaded / total * 100 if total else 0
                    speed = downloaded / (time.time() - start) / 1e6
                    print(f"  [{SAMPLE}] R{read_num}: {downloaded/1e9:.2f}/{total/1e9:.2f} GB ({pct:.1f}%) @ {speed:.1f} MB/s", flush=True)
        SIZE_TOTAL += downloaded
        print(f"  [{SAMPLE}] R{read_num} downloaded ({downloaded/1e9:.2f} GB)", flush=True)
    except Exception as e:
        print(f"  [{SAMPLE}] R{read_num} download FAILED: {type(e).__name__}", flush=True)
        sys.exit(1)

# Step 2: Run salmon quant in WSL
quant_dir = f"{OUTDIR_WSL}/{SAMPLE}_{COND}_{STAGE}"
r1 = f"{TMPDIR_WSL}/{SAMPLE}_1.fastq.gz"
r2 = f"{TMPDIR_WSL}/{SAMPLE}_2.fastq.gz"

if os.path.exists(f"{OUTDIR_WIN}/{SAMPLE}_{COND}_{STAGE}/quant.sf"):
    print(f"  [{SAMPLE}] Already quantified, skipping", flush=True)
else:
    print(f"  [{SAMPLE}] Running salmon quant...", flush=True)
    cmd = f'wsl --distribution Ubuntu --exec bash -c "salmon quant -i {INDEX_WSL} -l A -1 {r1} -2 {r2} --validateMappings -p 4 -o {quant_dir}"'
    rc = os.system(cmd)
    if rc != 0:
        print(f"  [{SAMPLE}] salmon quant FAILED (rc={rc})", flush=True)
        sys.exit(1)
    print(f"  [{SAMPLE}] salmon quant done", flush=True)

# Step 3: Clean up FASTQ files
for read_num in ["1", "2"]:
    p = f"{TMPDIR_WIN}/{SAMPLE}_{read_num}.fastq.gz"
    if os.path.exists(p):
        os.remove(p)
        print(f"  [{SAMPLE}] cleaned up R{read_num}", flush=True)

elapsed = time.time() - start
speed = SIZE_TOTAL / elapsed / 1e6 if elapsed > 0 else 0
print(f"  [{SAMPLE}] DONE in {elapsed/60:.1f} min ({speed:.1f} MB/s avg)", flush=True)
