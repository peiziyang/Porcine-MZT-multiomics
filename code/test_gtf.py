import gzip, time
f='E:/Workbuddy/2026-07-27-11-58-27/data/raw/pig_genes.gtf.gz'
start=time.time()
n=0; genes=0
with gzip.open(f,'rt') as fh:
    for line in fh:
        if line.startswith('#'): continue
        parts=line.strip().split('\t')
        if len(parts)>=9 and parts[2]=='gene':
            genes+=1
        n+=1
        if n%50000==0:
            print(f'  {n} lines, {genes} genes ({time.time()-start:.0f}s)', flush=True)
        if n>200000:
            break
print(f'Scanned {n} lines, {genes} gene annotations in {time.time()-start:.1f}s', flush=True)
