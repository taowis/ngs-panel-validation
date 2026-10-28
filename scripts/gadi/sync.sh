#!/bin/bash
# Copy report 05 scripts and inputs to Gadi, and bring the results back.
#   scripts/gadi/sync.sh push     scripts, capture BED, HG003 capture BAM, reference .fai
#   scripts/gadi/sync.sh pull     VCFs, logs and run records into work_dv/gadi/
# The BAM is built locally first:
#   samtools view -M -L work_dv/capture.HG003.bed -b -o work_dv/HG003.capture.bam results_v2_HG003/bam/HG003.md.bam
set -euo pipefail
cd "$(dirname "$0")/../.."
HOST=ktl590@gadi.nci.org.au
BASE=/scratch/iz89/ktl590/ngs-panel-validation

case $1 in
push)
    ssh -o BatchMode=yes $HOST "mkdir -p $BASE/scripts/gadi $BASE/inputs $BASE/ref"
    rsync -av scripts/gadi/ $HOST:$BASE/scripts/gadi/
    rsync -av --partial work_dv/capture.HG003.bed work_dv/HG003.capture.bam work_dv/HG003.capture.bam.bai $HOST:$BASE/inputs/
    rsync -av data/reference/hs37d5.fa.fai $HOST:$BASE/ref/
    ;;
pull)
    mkdir -p work_dv/gadi
    rsync -av --include='*/' --include='*.vcf.gz' --include='*.vcf.gz.tbi' --include='*.txt' --include='*.log' \
        --include='*.md5' --include='*.o*' --exclude='*' $HOST:$BASE/results $HOST:$BASE/logs "$HOST:$BASE/dv*.o*" work_dv/gadi/
    ;;
esac
