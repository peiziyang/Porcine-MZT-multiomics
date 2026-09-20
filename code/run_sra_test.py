"""Run v3 SRA test logging to file."""
import os, sys
cmd = f'C:/Users/peiya/.workbuddy/binaries/python/envs/default/Scripts/python.exe E:/Workbuddy/2026-07-27-11-58-27/scripts/sra_salmon_pipeline_v3.py SRR13435642 IVF 1cell > E:/Workbuddy/2026-07-27-11-58-27/data/raw/sra_download_test.log 2>&1'
print(f"Starting: {cmd[:80]}...")
rc = os.system(cmd)
print(f"Exit: {rc}")
