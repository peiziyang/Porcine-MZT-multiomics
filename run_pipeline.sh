#!/bin/bash
# Full pipeline: convert the 2 re-downloaded .sra -> fastq.gz (only the missing read each),
# then re-run salmon for the 2 failed PA 2cell samples, then finalize.
# TMP on E: (not C:) to avoid slow WSL rootfs VHDX and C: fill risk.
# ── CONFIG: edit these to match your environment (or export env vars) ──
BIN="${SRATOOLKIT_BIN:-/REPLACE_WITH_PATH/to/sratoolkit/bin}"   # dir containing fasterq-dump
SRADIR="${SRA_DIR:-/REPLACE_WITH_PATH/to/data/raw/sra_tmp}"     # where .sra files live
TMP="${TMP_DIR:-/REPLACE_WITH_PATH/to/tmp}"                     # temp space
LOGDIR="${LOG_DIR:-/REPLACE_WITH_PATH/to/data/processed/salmon_out}"
OUTDIR=$LOGDIR
PIPELOG=$LOGDIR/run_pipeline.log
exec > "$PIPELOG" 2>&1
echo "=== PIPELINE START $(date) ==="
mkdir -p "$TMP"

echo "--- convert 672 (keep _2) ---"
"$BIN/fasterq-dump" "$SRADIR/SRR13435672" -O "$TMP" --split-files -e 4 -p --temp "$TMP" > "$LOGDIR/conv672.log" 2>&1
echo "fq672 exit=$?"
gzip "$TMP/SRR13435672_2.fastq"
mv "$TMP/SRR13435672_2.fastq.gz" "$SRADIR/"
echo "672_2 moved"

echo "--- convert 673 (keep _1) ---"
"$BIN/fasterq-dump" "$SRADIR/SRR13435673" -O "$TMP" --split-files -e 4 -p --temp "$TMP" > "$LOGDIR/conv673.log" 2>&1
echo "fq673 exit=$?"
gzip "$TMP/SRR13435673_1.fastq"
mv "$TMP/SRR13435673_1.fastq.gz" "$SRADIR/"
echo "673_1 moved"

rm -rf "$TMP"
echo "--- integrity ---"
gzip -t "$SRADIR/SRR13435672_2.fastq.gz" && echo 672_2_OK || echo 672_2_BAD
gzip -t "$SRADIR/SRR13435673_1.fastq.gz" && echo 673_1_OK || echo 673_1_BAD
rm -f "$SRADIR/SRR13435672_2.fastq.gz.corrupt" "$SRADIR/SRR13435673_1.fastq.gz.corrupt"

echo "--- salmon rerun (2 samples) ---"
bash "$(dirname "$0")/rerun_failed.sh"

N=$(find "$OUTDIR" -name quant.sf | wc -l)
echo "quant.sf count = $N"
if [ "$N" -ge 42 ]; then
  echo "All done" >> "$OUTDIR/batch_wsl2.log"
  echo "=== PIPELINE COMPLETE 42/42 $(date) ==="
else
  echo "WARNING: only $N/42 quant.sf produced"
fi
