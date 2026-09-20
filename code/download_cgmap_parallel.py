"""Multi-threaded CGmap download with parallel HTTP range requests."""
import os
import sys
import time
import urllib.request
import threading

URL = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE235nnn/GSE235731/suppl/GSE235731_RAW.tar'
DEST = 'E:/Workbuddy/2026-07-27-11-58-27/data/raw/cgmap/GSE235731_RAW.tar'
N_PARALLEL = 4
CHUNK_SIZE = 8 * 1024 * 1024  # 8 MB

# Get total file size
req = urllib.request.Request(URL, headers={'User-Agent': 'Mozilla/5.0'})
resp = urllib.request.urlopen(req, timeout=30)
total = int(resp.headers.get('content-length', 0))
print(f'Total: {total/1e9:.1f} GB, using {N_PARALLEL} parallel threads')

# Check existing file
downloaded = 0
if os.path.exists(DEST + '.tmp'):
    downloaded = os.path.getsize(DEST + '.tmp')
    print(f'Resuming from {downloaded/1e9:.2f} GB')

# Calculate chunks
chunk_starts = list(range(downloaded, total, CHUNK_SIZE))
chunk_ends = [min(s + CHUNK_SIZE, total) for s in chunk_starts]

# Track progress
progress_lock = threading.Lock()
progress = [downloaded] * N_PARALLEL
last_log = time.time()
start_time = time.time()

def download_chunk(idx, start, end):
    global progress
    pos = start
    tmp = DEST + f'.part{idx}'
    try:
        while pos < end:
            req = urllib.request.Request(URL, headers={
                'User-Agent': 'Mozilla/5.0',
                'Range': f'bytes={pos}-{end-1}'
            })
            resp = urllib.request.urlopen(req, timeout=60)
            with open(tmp, 'ab') as f:
                while pos < end:
                    chunk = resp.read(1024 * 1024)  # 1 MB
                    if not chunk:
                        break
                    f.write(chunk)
                    pos += len(chunk)
                    with progress_lock:
                        progress[idx] = pos - start + (idx * (CHUNK_SIZE if idx else 0))
    except Exception as e:
        print(f'Thread {idx}: error {e}', flush=True)

# Check if server supports range requests
test_req = urllib.request.Request(URL, headers={
    'User-Agent': 'Mozilla/5.0',
    'Range': f'bytes={downloaded}-{downloaded + CHUNK_SIZE - 1}'
})
try:
    test_resp = urllib.request.urlopen(test_req, timeout=15)
    if test_resp.status == 206 or test_resp.status == 200:
        print('Server supports range requests - using parallel download')
    else:
        print(f'Unexpected status: {test_resp.status}')
except Exception as e:
    print(f'Range request failed: {e} - falling back to single thread')
    # Fall back to single thread download
    os.system(f'python {sys.argv[0] if len(sys.argv) > 1 else "download_cgmap_resume.py"}')

# Launch threads
threads = []
for i, (s, e) in enumerate(zip(chunk_starts, chunk_ends)):
    if s >= total:
        break
    t = threading.Thread(target=download_chunk, args=(i, s, e))
    t.start()
    threads.append(t)

# Wait for completion with progress logging
for t in threads:
    t.join()

# Combine parts
with open(DEST, 'wb') as out:
    for i in range(len(chunk_starts)):
        part = DEST + f'.part{i}'
        if os.path.exists(part):
            with open(part, 'rb') as f:
                while True:
                    chunk = f.read(8 * 1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
            os.remove(part)

elapsed = time.time() - start_time
print(f'Download complete in {elapsed/60:.1f} min ({os.path.getsize(DEST)/1e9:.2f} GB)')