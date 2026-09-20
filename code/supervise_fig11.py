"""
supervise_fig11.py — self-healing runner for the Fig 11 GO enrichment pipeline.

Why this exists:
  A script can "succeed" (exit 0) yet produce an EMPTY/WRONG result (e.g. 0 GO
  terms because a dict was keyed the wrong way). Plain background runs only alert
  on crashes, so a silent bad result only surfaces when the user asks next time.
  This supervisor closes that gap:

  1. runs m15_fig11_go_enrich.py via subprocess
  2. VALIDATES outputs (CSV non-empty, GO terms tested > 0, PNG > 1 KB)
  3. on failure: reads the log tail, matches a known error signature, applies a
     source-level patch to the script, and retries (up to MAX_TRIES)
  4. on success: runs build_v7_fig11.py, regenerates v7.docx, validates images
  5. writes SUPERVISE_REPORT.txt and a clear final status

Run it in background; it notifies on completion. No manual "is it done?" needed.
"""
import os, re, subprocess, sys, time, zipfile

ROOT   = 'E:/Workbuddy/2026-07-27-11-58-27'
SCRIPT = os.path.join(ROOT, 'scripts', 'm15_fig11_go_enrich.py')
BUILD  = os.path.join(ROOT, 'scripts', 'build_v7_fig11.py')
DOCX   = os.path.join(ROOT, 'scripts', 'md_to_docx.py')
CSV    = os.path.join(ROOT, 'data/processed/deseq2', 'm15_fig11_go_enrich.csv')
PNG    = os.path.join(ROOT, 'data/processed/deseq2', 'm15_fig11_go_enrich.png')
LOG    = os.path.join(ROOT, 'data/processed/deseq2', 'm15_fig11.log')
V7MD   = os.path.join(ROOT, 'manuscript', 'manuscript_full_with_figures_v7.md')
V7DOCX = os.path.join(ROOT, 'manuscript', 'manuscript_full_with_figures_v7.docx')
REPORT = os.path.join(ROOT, 'data/processed/deseq2', 'SUPERVISE_REPORT.txt')
PY     = 'C:/Users/peiya/.workbuddy/binaries/python/envs/default/Scripts/python.exe'
MAX_TRIES = 5

def log(msg):
    print(msg)
    with open(REPORT, 'a', encoding='utf-8') as f:
        f.write(msg + '\n')

def read_log_tail(n=50):
    try:
        lines = open(LOG, encoding='utf-8', errors='ignore').read().splitlines()
        return '\n'.join(lines[-n:])
    except Exception:
        return ''

def count_sig_terms():
    """Number of GO terms with BH q < 0.05 (significant). -1 if unreadable."""
    try:
        import pandas as pd
        df = pd.read_csv(CSV)
        for c in ['q', 'p']:
            if c in df.columns:
                v = pd.to_numeric(df[c], errors='coerce')
                return int((v < 0.05).sum())
    except Exception:
        return -1
    return -1

def validate():
    """Return (status_ok, detail).
    status_ok True  -> has significant terms; safe to merge.
    status_ok False + 'REVIEW' in detail -> ran OK but 0 significant terms
        (scientific null / low power). NOT a code bug -> do NOT auto-retry/build.
    status_ok False (no REVIEW) -> technical failure -> diagnose/patch/retry.
    """
    if not os.path.exists(CSV):
        return False, 'CSV missing'
    try:
        import pandas as pd
        df = pd.read_csv(CSV)
        if len(df) == 0:
            return False, 'CSV empty (0 rows)'
    except Exception as e:
        return False, f'CSV unreadable: {e}'
    if not os.path.exists(PNG) or os.path.getsize(PNG) < 1024:
        return False, 'PNG missing/too small'
    tail = read_log_tail()
    m = re.search(r'GO terms tested: (\d+)', tail)
    if m and int(m.group(1)) == 0:
        return False, 'GO terms tested = 0 (empty enrichment)'
    nsig = count_sig_terms()
    if nsig == 0:
        return False, ('REVIEW: 0 significant GO terms after BH (q<0.05); '
                      'likely true null or low power — needs human decision, NOT auto-retried')
    return True, f'ok ({nsig} significant terms)'

# --- known source-level patches (signature -> patch fn) -----------------------
def patch_go2genes():
    src = open(SCRIPT, encoding='utf-8').read()
    if 'go2genes' in src and 'for go, term_genes in go2genes.items()' in src:
        return False  # already applied
    src = src.replace(
        "print(f'    GO terms observed: {len(go_info)}')",
        "print(f'    GO terms observed: {len(go_info)}')\n\n"
        "# INVERT: go -> gene set for enrichment (entrez2go above is gene->GO)\n"
        "go2genes = {}\n"
        "for _gid, _gos in entrez2go.items():\n"
        "    for _go in _gos:\n"
        "        go2genes.setdefault(_go, set()).add(_gid)")
    src = src.replace('for go, term_genes in entrez2go.items():',
                      'for go, term_genes in go2genes.items():')
    open(SCRIPT, 'w', encoding='utf-8').write(src)
    return True

def patch_tonumeric():
    src = open(SCRIPT, encoding='utf-8').read()
    if "pd.to_numeric(res['p']" in src:
        return False
    src = src.replace(
        "res = res[np.isfinite(res['p'])].copy()",
        "res['p'] = pd.to_numeric(res['p'], errors='coerce')\n"
        "res = res[np.isfinite(res['p'])].copy()")
    open(SCRIPT, 'w', encoding='utf-8').write(src)
    return True

def patch_subset():
    src = open(SCRIPT, encoding='utf-8').read()
    if 'deg_with_go = deg_with_go & bg_with_go' in src:
        return False
    src = src.replace(
        "deg_with_go = deg_entrez & set(entrez2go.keys())",
        "deg_with_go = deg_entrez & set(entrez2go.keys())\n"
        "deg_with_go = deg_with_go & bg_with_go  # ensure N<=M for hypergeom")
    open(SCRIPT, 'w', encoding='utf-8').write(src)
    return True

PATCHES = [
    (r'GO terms tested: 0',               patch_go2genes, 'invert entrez2go -> go2genes'),
    (r"ufunc 'isfinite' not supported",   patch_tonumeric, 'to_numeric before isfinite'),
    (r'must contain only real numbers',   patch_tonumeric, 'to_numeric before FDR'),
    (r'The count matrix should only',     None,           'count matrix must be integer'),
]

# --- main loop ---------------------------------------------------------------
open(REPORT, 'w', encoding='utf-8').write('SUPERVISE FIG11 — start\n')
ok = False
review_mode = False
for attempt in range(1, MAX_TRIES + 1):
    log(f'[try {attempt}] running enrichment...')
    t0 = time.time()
    r = subprocess.run([PY, SCRIPT], capture_output=True, text=True)
    with open(LOG, 'w', encoding='utf-8') as f:
        f.write(r.stdout + '\n' + r.stderr)
    log(f'[try {attempt}] exit={r.returncode} in {time.time()-t0:.1f}s')
    ok, detail = validate()
    if ok:
        log(f'[try {attempt}] VALIDATION OK: {detail}')
        break
    if 'REVIEW' in detail:
        log(f'[try {attempt}] VALIDATION REVIEW: {detail}')
        review_mode = True
        break
    log(f'[try {attempt}] VALIDATION FAIL: {detail}')
    # diagnose: match signature in log tail or stderr
    text = read_log_tail() + '\n' + r.stderr
    applied = False
    for sig, fn, desc in PATCHES:
        if re.search(sig, text):
            if fn is None:
                log(f'[try {attempt}] known issue but no auto-patch ({desc}); stopping')
                applied = False
                break
            if fn():
                log(f'[try {attempt}] applied patch: {desc}')
                applied = True
                break
    if not applied:
        log(f'[try {attempt}] no applicable auto-patch; stopping to avoid blind retry')
        break

if review_mode:
    log('FINAL: NEEDS_REVIEW — enrichment ran but 0 BH-significant GO terms. '
        'Treat as scientific null / low power, NOT a code bug. Do NOT auto-merge; '
        'manually revise v7 text (the raw-p figure is informative) then regenerate docx.')
    sys.exit(0)
if not ok:
    log('FINAL: enrichment failed after attempts — see m15_fig11.log and SUPERVISE_REPORT.txt')
    sys.exit(2)

# --- success path: merge into v7 + regenerate docx ---------------------------
log('[build] merging Fig 11 into v7.md ...')
rb = subprocess.run([PY, BUILD], capture_output=True, text=True)
log(rb.stdout.strip())
log(rb.stderr.strip())
v7 = open(V7MD, encoding='utf-8').read()
if '## Figure 11' not in v7:
    log('[build] ERROR: Fig 11 section not found in v7.md'); sys.exit(3)
log('[docx] regenerating v7.docx ...')
rd = subprocess.run([PY, DOCX, V7MD, V7DOCX], capture_output=True, text=True)
log(rd.stdout.strip()); log(rd.stderr.strip())
try:
    nimg = len([n for n in zipfile.ZipFile(V7DOCX).namelist() if n.startswith('word/media/')])
    log(f'[docx] embedded images = {nimg}')
except Exception as e:
    log(f'[docx] could not inspect docx: {e}')
log('FINAL: SUCCESS — Fig 11 enrichment + v7 integrated and docx regenerated.')
