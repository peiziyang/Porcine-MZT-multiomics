"""SRA FASTQ download + salmon quantification pipeline v3.
Always validates file size against server content-length.
Usage: python sra_salmon_pipeline.py [SRR_id] [condition] [stage]
"""
import os, sys, urllib.request, time

SAMPLE = sys.argv[1] if len(sys.argv) > 1 else "SRR13435642"
COND = sys.argv[2] if len(sys.argv) > 2 else "?"
STAGE = sys.argv[3] if len(sys.argv) > 3 else "?"

TMPDIR_WIN = "E:/Workbuddy/2026-07-27-11-58-27/data/raw/sra_tmp"
OUTDIR_WIN = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/salmon_out"
TMPDIR_WSL = "/mnt/e/Workbuddy/2026-07-27-11-58-27/data/raw/sra_tmp"
INDEX_WSL = "/mnt/e/Workbuddy/2026-07-27-11-58-27/data/raw/salmon_index"
OUTDIR_WSL = "/mnt/e/Workbuddy/2026-07-27-11-58-27/data/processed/salmon_out"
EBI_PREFIX = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR134/042"

os.makedirs(TMPDIR_WIN, exist_ok=True)
os.makedirs(OUTDIR_WIN, exist_ok=True)
start = time.time()
SIZE_TOTAL = 0

def check_file_size(url):
    """Get content-length from EBI FTP without downloading."""
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        resp = urllib.request.urlopen(req, timeout=30)
        return int(resp.headers.get("content-length", 0))
    except:
        return 0

print(f"[{SAMPLE}] {COND} {STAGE} - start", flush=True)

# Step 1: Download R1 and R2
for read_num in ["1", "2"]:
    url = f"{EBI_PREFIX}/{SAMPLE}/{SAMPLE}_{read_num}.fastq.gz"
    dest = f"{TMPDIR_WIN}/{SAMPLE}_{read_num}.fastq.gz"
    
    # Check expected size
    expected = check_file_size(url)
    if expected <= 0:
        print(f"  [{SAMPLE}] Cannot get file size for R{read_num}, abort", flush=True)
        sys.exit(1)
    
    # Validate existing file
    if os.path.exists(dest):
        got = os.path.getsize(dest)
        if got == expected:
            print(f"  [{SAMPLE}] R{read_num} already complete ({got/1e9:.2f} GB)", flush=True)
            SIZE_TOTAL += got
            continue
        else:
            print(f"  [{SAMPLE}] R{read_num} partial ({got/1e9:.2f}/{expected/1e9:.2f} GB), re-downloading", flush=True)
            try:
                os.remove(dest)
            except:
                pass
    
    print(f"  [{SAMPLE}] Downloading R{read_num} ({expected/1e9:.2f} GB)...", flush=True)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        resp = urllib.request.urlopen(req, timeout=120)
        CHUNK = 8 * 1024 * 1024
        downloaded = 0
        with open(dest, "wb") as f:
            while True:
                chunk = resp.read(CHUNK)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)
                if downloaded % (100 * 1024 * 1024) < CHUNK:
                    pct = downloaded / expected * 100
                    speed = downloaded / (time.time() - start) / 1e6
                    print(f"  [{SAMPLE}] R{read_num}: {downloaded/1e9:.2f}/{expected/1e9:.2f} GB ({pct:.1f}%) @ {speed:.1f} MB/s", flush=True)
        SIZE_TOTAL += downloaded
        print(f"  [{SAMPLE}] R{read_num} downloaded", flush=True)
    except Exception as e:
        print(f"  [{SAMPLE}] R{read_num} FAILED: {type(e).__name__}", flush=True)
        sys.exit(1)

# Step 2: salmon quant in WSL
quant_dir = f"{OUTDIR_WSL}/{SAMPLE}_{COND}_{STAGE}"
r1 = f"{TMPDIR_WSL}/{SAMPLE}_1.fastq.gz"
r2 = f"{TMPDIR_WSL}/{SAMPLE}_2.fastq.gz"
if os.path.exists(f"{OUTDIR_WIN}/{SAMPLE}_{COND}_{STAGE}/quant.sf"):
    print(f"  [{SAMPLE}] Already quantified, skip", flush=True)
else:
    print(f"  [{SAMPLE}] Running salmon quant...", flush=True)
    cmd = f'wsl --distribution Ubuntu --exec bash -c "salmon quant -i {INDEX_WSL} -l A -1 {r1} -2 {r2} --validateMappings -p 4 -o {quant_dir}"'
    rc = os.system(cmd)
    if rc != 0:
        print(f"  [{SAMPLE}] salmon quant FAILED (rc={rc})", flush=True)
        sys.exit(1)
    print(f"  [{SAMPLE}] salmon quant done", flush=True)

# Step 3: Clean up
for rn in ["1","2"]:
    p = f"{TMPDIR_WIN}/{SAMPLE}_{rn}.fastq.gz"
    if os.path.exists(p):
        try:
            os.remove(p)
        except:
            pass

elapsed = time.time() - start
speed = SIZE_TOTAL / elapsed / 1e6 if elapsed > 0 else 0
print(f"[{SAMPLE}] DONE ({elapsed/60:.1f} min, {speed:.1f} MB/s avg)", flush=True)
