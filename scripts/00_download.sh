#!/usr/bin/env bash
# Fetch public inputs. URLs are pinned to fixed releases; checksums are in
# data/MANIFEST.tsv and checked after download.
#   scripts/00_download.sh panel   annotation, capture BED, truth set (~200 MB)
#   scripts/00_download.sh reference  hs37d5 FASTA (~900 MB) and indexes
#   scripts/00_download.sh reads   HG001 exome FASTQ (~30 GB)
set -euo pipefail
cd "$(dirname "$0")/.."

giab=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab
garvan=$giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome
truth=$giab/release/NA12878_HG001/NISTv4.2.1/GRCh37
# GIAB's AWS mirror holds the same files (same sizes) and is about 50x faster
# than ftp-trace from Sydney, so the large FASTQ come from there.
garvan_s3=https://giab.s3.amazonaws.com/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome

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
        fetch $garvan_s3/${lib}_${lane}_${r}_001.fastq.gz data/raw
      done
    done
  done
}

reference() {
  # 1000 Genomes AWS mirror: the original 2011 upload. EBI serves a 2017
  # recompression (different .gz size) with an identical .fai, and is far slower.
  local ref=https://1000genomes.s3.amazonaws.com/technical/reference/phase2_reference_assembly_sequence
  fetch $ref/hs37d5.fa.gz data/reference
  # derived: rebuilt from the download above
  # the 2011 file ends in an index block that gzip reports as trailing garbage
  # (exit 2); accept that only if the output has the size the .fai implies
  if [ ! -s data/reference/hs37d5.fa ]; then
    gzip -dc data/reference/hs37d5.fa.gz > data/reference/hs37d5.fa || [ $? -eq 2 ]
    [ "$(wc -c < data/reference/hs37d5.fa | tr -d ' ')" = 3189750467 ]
  fi
  [ -s data/reference/hs37d5.fa.fai ] || samtools faidx data/reference/hs37d5.fa
  [ -s data/reference/hs37d5.dict ] || gatk CreateSequenceDictionary -R data/reference/hs37d5.fa
  [ -s data/reference/hs37d5.fa.bwt ] || bwa index data/reference/hs37d5.fa
  [ -d data/reference/hs37d5.sdf ] || rtg format -o data/reference/hs37d5.sdf data/reference/hs37d5.fa
}

"${1:-panel}"
python3 scripts/check_manifest.py
