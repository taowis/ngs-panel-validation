#!/usr/bin/env bash
# Index the test reference (the same indexes scripts/00_download.sh builds for hs37d5).
set -euo pipefail
cd "$(dirname "$0")/data/ref"
[ -s hs37d5.fa.fai ] || samtools faidx hs37d5.fa
[ -s hs37d5.dict ] || gatk CreateSequenceDictionary -R hs37d5.fa -O hs37d5.dict --QUIET true 2> /dev/null
[ -s hs37d5.fa.bwt ] || bwa index hs37d5.fa 2> /dev/null
[ -d hs37d5.sdf ] || rtg format -o hs37d5.sdf hs37d5.fa > /dev/null
