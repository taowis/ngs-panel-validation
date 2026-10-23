#!/usr/bin/env bash
# Run the pipeline for each genome in turn (report 03), one output folder each.
# Builds the genome's panel evaluable region first. Start it detached:
#   nohup scripts/run_genomes.sh HG002 HG003 HG004 HG005 > results_v2_genomes.log 2>&1 < /dev/null &
# -resume picks up the most recent Nextflow session, so after an interruption
# rerun the interrupted genome first, then the ones not yet run.
set -euo pipefail
cd "$(dirname "$0")/.."

for g in "$@"; do
  p=config/params/$g.yaml
  val() { awk -v k="$1:" '$1 == k { print $2 }' "$p"; }
  out=$(val outdir)
  mkdir -p "$out/pipeline_info"
  python3 scripts/01_build_panel.py --tag "$g" --capture "$(val capture)" \
    --truth-bed "$(val truth_bed)" --truth-vcf "$(val truth_vcf)" > /dev/null
  echo "$(date '+%F %T') $g start"
  { git rev-parse HEAD; caffeinate -i nextflow run workflow/main.nf -profile local -resume -ansi-log false \
      -params-file "$p"; } > "$out/pipeline_info/run.log" 2>&1
  echo "$(date '+%F %T') $g done"
done
