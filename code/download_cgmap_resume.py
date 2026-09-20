"""Restart CGmap download with HTTP range request for resume."""
import urllib.request, sys, os, time

url = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE235nnn/GSE235731/suppl/GSE235731_RAW.tar'
dest = 'E:/Workbuddy/2026-07-27-11-58-27/data/raw/cgmap/GSE235731_RAW.tar'
tmp = dest + '.tmp'

# Get total file size
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
resp = urllib.request.urlopen(req, timeout=30)
total = int(resp.headers.get('content-length', 0))
print(f'Total size: {total/1e9:.1f} GB', flush=True)

# Check existing partial file
if os.path.exists(tmp):
    downloaded = os.path.getsize(tmp)
    print(f'Resuming from {downloaded/1e9:.2f} GB', flush=True)
else:
    downloaded = 0
    print(f'Starting fresh download', flush=True)

CHUNK = 8 * 1024 * 1024  # 8 MB
last_log = time.time()
start_time = time.time()

while downloaded < total:
    try:
        req = urllib.request.Request(url, headers={
            'User-Agent': 'Mozilla/5.0',
            'Range': f'bytes={downloaded}-'
        })
        resp = urllib.request.urlopen(req, timeout=120)
    except Exception as e:
        print(f'Network error at {downloaded/1e9:.2f} GB: {type(e).__name__}, retry in 10s', flush=True)
        time.sleep(10)
        continue

    mode = 'ab' if downloaded > 0 else 'wb'
    with open(tmp, mode) as f:
        while True:
            chunk = resp.read(CHUNK)
            if not chunk:
                break
            f.write(chunk)
            downloaded += len(chunk)
            if time.time() - last_log > 60:
                pct = downloaded / total * 100
                elapsed = time.time() - start_time
                speed_mbps = downloaded / 1024 / 1024 / max(elapsed, 1)
                eta_min = (total - downloaded) / 1024 / 1024 / max(speed_mbps, 0.01)
                print(f'  [{time.strftime("%H:%M:%S")}] {downloaded/1e9:.2f}/{total/1e9:.1f} GB ({pct:.1f}%) | {speed_mbps:.1f} MB/s | ETA {eta_min/60:.1f}h', flush=True)
                last_log = time.time()
            if downloaded >= total:
                break

if os.path.exists(tmp):
    os.rename(tmp, dest)
    print(f'Done: {dest}', flush=True)