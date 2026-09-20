#!/bin/bash
# Re-run salmon for the 2 corrupted PA 2cell samples (SRR13435672, SRR13435673).
# Run inside WSL: wsl -d Ubuntu-22.04 bash ./rerun_failed.sh
# Prerequisite: user must re-download the 2 corrupted fastq.gz into sra_tmp:
#   SRR13435672_2.fastq.gz  (R2 of 672)
#   SRR13435673_1.fastq.gz  (R1 of 673)
SRADIR=${SRA_DIR:-/REPLACE_WITH_PATH/to/data/raw/sra_tmp}
INDEX=/home/salmon_index
OUTDIR=${OUT_DIR:-/REPLACE_WITH_PATH/to/data/processed/salmon_out}
LOGFILE=$OUTDIR/rerun_failed.log
exec > "$LOGFILE" 2>&1
echo "Rerun failed (2 samples): $(date)"
SAMPLES=(
  "SRR13435672 PA 2cell"
  "SRR13435673 PA 2cell"
)
for entry in "${SAMPLES[@]}"; do
  srr=$(echo $entry | cut -d' ' -f1)
  cond=$(echo $entry | cut -d' ' -f2)
  stage=$(echo $entry | cut -d' ' -f3)
  odir="$OUTDIR/${srr}_${cond}_${stage}"
  # sanity: both fastq must exist
  if [ ! -f "$SRADIR/${srr}_1.fastq.gz" ] || [ ! -f "$SRADIR/${srr}_2.fastq.gz" ]; then
    echo "$srr: MISSING fastq, abort this sample"
    continue
  fi
  rm -rf "$odir"
  echo "$srr ($cond $stage): start $(date)"
  salmon quant -i $INDEX -l A -1 $SRADIR/${srr}_1.fastq.gz -2 $SRADIR/${srr}_2.fastq.gz --validateMappings -p 2 -o "$odir"
  if [ -f "$odir/quant.sf" ]; then echo "$srr: done"; else echo "$srr: FAIL"; fi
done
echo "Rerun finished: $(date)"
