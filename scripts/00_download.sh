#!/usr/bin/env bash
# Fetch public inputs. URLs are pinned to fixed releases; checksums are in
# data/MANIFEST.tsv and checked after download.
#   scripts/00_download.sh panel   annotation, capture BED, truth set (~200 MB)
#   scripts/00_download.sh reads   HG001 exome FASTQ (~30 GB, run on HPC)
set -euo pipefail
cd "$(dirname "$0")/.."

giab=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab
garvan=$giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome
truth=$giab/release/NA12878_HG001/NISTv4.2.1/GRCh37

fetch() {  # url dest_dir
  mkdir -p "$2"
  local f="$2/$(basename "$1")"
  if [ ! -s "$f" ]; then
    curl -sSfL --retry 3 -o "$f.part" "$1"
    mv "$f.part" "$f"
  fi
}

panelapp() {
  # PanelApp Australia API, pinned panel versions. Responses are committed in
  # data/panels/panelapp/ because the API can change; rerun only to update.
  mkdir -p data/panels/panelapp
  for p in 4375:2.0 4374:2.0 4371:2.1; do
    curl -sSf "https://panelapp-aus.org/api/v1/panels/${p%%:*}/?version=${p##*:}" \
      | python3 -m json.tool --sort-keys > "data/panels/panelapp/panel_${p%%:*}_v${p##*:}.json"
  done
}

panel() {
  fetch https://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz data/reference
  fetch $garvan/nexterarapidcapture_expandedexome_targetedregions.bed.gz data/raw
  for f in HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi HG001_GRCh37_1_22_v4.2.1_benchmark.bed; do
    fetch $truth/$f data/truth
  done
}

reads() {
  for lib in NIST7035_TAAGGCGA NIST7086_CGTACTAG; do
    for lane in L001 L002; do
      for r in R1 R2; do
        fetch $garvan/${lib}_${lane}_${r}_001.fastq.gz data/raw
      done
    done
  done
}

"${1:-panel}"
python3 scripts/check_manifest.py
