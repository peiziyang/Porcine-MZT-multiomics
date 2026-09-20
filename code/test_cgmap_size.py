"""Estimate CGmap file size and CpG count."""
import gzip, time, os
f='E:/Workbuddy/2026-07-27-11-58-27/data/raw/cgmap/GSM7508586_AF1_12.CGmap.gz'
start=time.time()
total=0; cpg=0
with gzip.open(f,'rt') as fh:
    for i,line in enumerate(fh):
        total+=1
        if line.split('\t')[3]=='CG':
            cpg+=1
        if i>=500000:
            break
size=os.path.getsize(f)
ratio=size/(i+1) if i>0 else 1
est_lines=size/ratio
print(f'Sample scan: {i+1} lines ({cpg} CpG) in {time.time()-start:.0f}s')
print(f'Estimated total lines: {est_lines/1e6:.0f}M, CpG: {cpg/(i+1)*100:.1f}%')
print(f'Estimated CpG sites per file: {est_lines*cpg/(i+1)/1e6:.0f}M')
