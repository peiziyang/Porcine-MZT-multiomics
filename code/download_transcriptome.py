"""Download pig transcriptome FASTA from Ensembl for salmon index building."""
import urllib.request, sys, os

URL = "https://ftp.ensembl.org/pub/release-108/fasta/sus_scrofa/cdna/Sus_scrofa.Sscrofa11.1.cdna.all.fa.gz"
DEST = "E:/Workbuddy/2026-07-27-11-58-27/data/raw/pig_transcriptome.fa.gz"

req = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0"})
resp = urllib.request.urlopen(req, timeout=120)
total = int(resp.headers.get("content-length", 0))
print(f"Transcriptome size: {total/1e6:.1f} MB")

CHUNK = 4 * 1024 * 1024
downloaded = 0
with open(DEST, "wb") as f:
    while True:
        chunk = resp.read(CHUNK)
        if not chunk:
            break
        f.write(chunk)
        downloaded += len(chunk)
        pct = downloaded / total * 100 if total else 0
        if downloaded % (10 * 1024 * 1024) < CHUNK:
            print(f"  {downloaded/1e6:.1f}/{total/1e6:.1f} MB ({pct:.0f}%)")

print(f"Download complete: {DEST}")
