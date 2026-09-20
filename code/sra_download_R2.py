"""Generate correct SRA download lists and start R2 batch download."""
import os

SRRS = [
    ("SRR13435642","IVF","1cell"),("SRR13435641","IVF","1cell"),("SRR13435640","IVF","1cell"),
    ("SRR13435660","IVF","2cell"),("SRR13435659","IVF","2cell"),("SRR13435658","IVF","2cell"),
    ("SRR13435657","IVF","4cell"),("SRR13435656","IVF","4cell"),("SRR13435655","IVF","4cell"),
    ("SRR13435654","IVF","8cell"),("SRR13435653","IVF","8cell"),("SRR13435652","IVF","8cell"),
    ("SRR13435651","IVF","8cell"),("SRR13435650","IVF","8cell"),("SRR13435649","IVF","8cell"),
    ("SRR13435648","IVF","4cell"),("SRR13435647","IVF","2cell"),("SRR13435646","IVF","1cell"),
    ("SRR13435645","IVF","2cell"),("SRR13435644","IVF","4cell"),("SRR13435643","IVF","1cell"),
    ("SRR13435681","PA","1cell"),("SRR13435680","PA","1cell"),("SRR13435679","PA","1cell"),
    ("SRR13435678","PA","1cell"),("SRR13435677","PA","1cell"),("SRR13435676","PA","1cell"),
    ("SRR13435675","PA","2cell"),("SRR13435674","PA","2cell"),("SRR13435673","PA","2cell"),
    ("SRR13435672","PA","2cell"),("SRR13435671","PA","2cell"),("SRR13435670","PA","2cell"),
    ("SRR13435669","PA","4cell"),("SRR13435668","PA","4cell"),("SRR13435667","PA","4cell"),
    ("SRR13435666","PA","4cell"),("SRR13435665","PA","4cell"),("SRR13435664","PA","4cell"),
    ("SRR13435663","PA","8cell"),("SRR13435662","PA","8cell"),("SRR13435661","PA","8cell"),
]

TMPDIR = "/mnt/e/Workbuddy/2026-07-27-11-58-27/data/raw/sra_tmp"

# Start R2 downloads on WSL
import subprocess, time

for srr, cond, stage in SRRS:
    # EBI directory = last 2 digits of SRR, zero-padded
    last2 = srr[-2:]
    dir3 = last2.zfill(3)  # e.g., "42" → "042"
    
    # R2 URL
    url = f"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR134/{dir3}/{srr}/{srr}_2.fastq.gz"
    dest = f"{TMPDIR}/{srr}_2.fastq.gz"
    
    # Skip if exists
    if os.path.exists(dest) and os.path.getsize(dest) > 1000000:
        print(f"[SKIP] {srr} R2 exists ({os.path.getsize(dest)/1e6:.0f} MB)", flush=True)
        continue
    
    cmd = f'wsl --distribution Ubuntu --exec bash -c "wget -q -O {dest} {url}"'
    print(f"[DL] {srr} R2 starting...", flush=True)
    subprocess.Popen(cmd, shell=True)
    time.sleep(1)  # stagger starts

print("All R2 downloads initiated!")
