"""Test EBI FTP URL patterns for different SRR accessions."""
import urllib.request

test_urls = [
    ("SRR13435642 (042)", "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR134/042/SRR13435642/SRR13435642_1.fastq.gz"),
    ("SRR13435641 (041)", "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR134/041/SRR13435641/SRR13435641_1.fastq.gz"),
    ("SRR13435681 (081)", "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR134/081/SRR13435681/SRR13435681_1.fastq.gz"),
    ("SRR13435660 (060)", "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR134/060/SRR13435660/SRR13435660_1.fastq.gz"),
]

for name, url in test_urls:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        resp = urllib.request.urlopen(req, timeout=15)
        size = int(resp.headers.get("content-length", 0))
        print(f"✅ {name}: {size/1e9:.2f} GB")
    except Exception as e:
        print(f"❌ {name}: {type(e).__name__}")
