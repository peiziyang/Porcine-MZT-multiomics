#!/bin/bash
# Salmon batch for all 42 GSE164812 samples (fresh WSL)
SRADIR=${SRA_DIR:-/REPLACE_WITH_PATH/to/data/raw/sra_tmp}
INDEX=/home/salmon_index
OUTDIR=${OUT_DIR:-/REPLACE_WITH_PATH/to/data/processed/salmon_out}
LOGFILE=$OUTDIR/batch_wsl2.log
mkdir -p $OUTDIR
exec > $LOGFILE 2>&1
echo "WSL2 salmon batch: $(date)"

SAMPLES=(
"SRR13435642 IVF 1cell" "SRR13435641 IVF 1cell" "SRR13435640 IVF 1cell" "SRR13435646 IVF 1cell" "SRR13435643 IVF 1cell"
"SRR13435660 IVF 2cell" "SRR13435659 IVF 2cell" "SRR13435658 IVF 2cell" "SRR13435647 IVF 2cell" "SRR13435645 IVF 2cell"
"SRR13435657 IVF 4cell" "SRR13435656 IVF 4cell" "SRR13435655 IVF 4cell" "SRR13435648 IVF 4cell" "SRR13435644 IVF 4cell"
"SRR13435654 IVF 8cell" "SRR13435653 IVF 8cell" "SRR13435652 IVF 8cell" "SRR13435651 IVF 8cell" "SRR13435650 IVF 8cell" "SRR13435649 IVF 8cell"
"SRR13435681 PA 1cell" "SRR13435680 PA 1cell" "SRR13435679 PA 1cell" "SRR13435678 PA 1cell" "SRR13435677 PA 1cell" "SRR13435676 PA 1cell"
"SRR13435675 PA 2cell" "SRR13435674 PA 2cell" "SRR13435673 PA 2cell" "SRR13435672 PA 2cell" "SRR13435671 PA 2cell" "SRR13435670 PA 2cell"
"SRR13435669 PA 4cell" "SRR13435668 PA 4cell" "SRR13435667 PA 4cell" "SRR13435666 PA 4cell" "SRR13435665 PA 4cell" "SRR13435664 PA 4cell"
"SRR13435663 PA 8cell" "SRR13435662 PA 8cell" "SRR13435661 PA 8cell"
)

run_one() {
    local entry="$1"
    local srr=$(echo $entry | cut -d' ' -f1)
    local cond=$(echo $entry | cut -d' ' -f2)
    local stage=$(echo $entry | cut -d' ' -f3)
    local st=$(date +%s)
    local odir="$OUTDIR/${srr}_${cond}_${stage}"
    if [ -f "$odir/quant.sf" ]; then echo "$srr: skip"; return; fi
    echo "$srr ($cond $stage): start $(date)"
    salmon quant -i $INDEX -l A -1 $SRADIR/${srr}_1.fastq.gz -2 $SRADIR/${srr}_2.fastq.gz --validateMappings -p 2 -o $odir 2>/dev/null
    local et=$(($(date +%s) - st))
    if [ -f "$odir/quant.sf" ]; then echo "$srr: done (${et}s)"; else echo "$srr: FAIL (${et}s)"; fi
}

N_PAR=4; idx=0; total=${#SAMPLES[@]}
while [ $idx -lt $total ]; do
    for ((i=0; i<N_PAR && idx+i<total; i++)); do
        run_one "${SAMPLES[$((idx+i))]}" &
    done
    idx=$((idx + N_PAR))
    wait
done
echo "All done: $(date)"
