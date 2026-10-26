#!/usr/bin/env bash
# Fetch public inputs. URLs are pinned to fixed releases; checksums are in
# data/MANIFEST.tsv and checked after download.
#   scripts/00_download.sh panel   annotation, capture BED, truth set (~200 MB)
#   scripts/00_download.sh reference  hs37d5 FASTA (~900 MB) and indexes
#   scripts/00_download.sh reads   HG001 exome FASTQ (~16 GB)
#   scripts/00_download.sh oslo    HG002-HG005 exome BAMs (~38 GB), FASTQ, truth sets
#   scripts/00_download.sh extra   ClinVar, GIAB stratifications, Oslo's own calls (~350 MB)
#   scripts/00_download.sh dbsnp   dbSNP 138 for the run QC known-variant rate (~1.5 GB)
set -euo pipefail
cd "$(dirname "$0")/.."

giab=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab
garvan=$giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome
truth=$giab/release/NA12878_HG001/NISTv4.2.1/GRCh37
# GIAB's AWS mirror holds the same files (same sizes) and is about 50x faster
# than ftp-trace from Sydney, so the large FASTQ come from there.
giab_s3=https://giab.s3.amazonaws.com
garvan_s3=$giab_s3/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome

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

oslo() {
  # GIAB HG002-HG005 exomes from Oslo University Hospital: Agilent SureSelect
  # Human All Exon V5, HiSeq 2500, one lane each. GIAB publishes only the aligned
  # BAMs, so the reads are taken back to FASTQ and realigned by the pipeline.
  local d=data/raw/oslo
  local a=$giab_s3/data/AshkenazimTrio/analysis/OsloUniversityHospital_Exome_GATK_jointVC_11242015
  fetch $a/wex_Agilent_SureSelect_v05_b37.baits.slop50.merged.list $d
  # derived: Picard interval list (1-based, closed) to BED
  if [ ! -s $d/sureselect_v5_b37.baits.slop50.bed.gz ]; then
    grep -v '^@' $d/wex_Agilent_SureSelect_v05_b37.baits.slop50.merged.list \
      | awk 'BEGIN { OFS = "\t" } { print $1, $2 - 1, $3 }' | gzip -n > $d/sureselect_v5_b37.baits.slop50.bed.gz
  fi
  local s g idx bam
  for s in AshkenazimTrio/HG002_NA24385_son:AGATGTAC AshkenazimTrio/HG003_NA24149_father:TCTTCACA \
           AshkenazimTrio/HG004_NA24143_mother:CCGAAGTA ChineseTrio/HG005_NA24631_son:CGCATACA; do
    g=$(basename ${s%%:*}); g=${g%%_*}; idx=${s##*:}
    bam=151002_7001448_0359_AC7F6GANXX_Sample_${g}-EEogPU_v02-KIT-Av5_${idx}_L008.posiSrt.markDup.bam
    fetch $giab_s3/data/${s%%:*}/OsloUniversityHospital_Exome/$bam $d &
  done
  wait
  for s in AshkenazimTrio/HG002_NA24385_son AshkenazimTrio/HG003_NA24149_father \
           AshkenazimTrio/HG004_NA24143_mother ChineseTrio/HG005_NA24631_son; do
    g=$(basename $s); g=${g%%_*}
    local t=$giab_s3/release/$s/NISTv4.2.1/GRCh37/${g}_GRCh37_1_22_v4.2.1_benchmark
    fetch $t.vcf.gz data/truth
    fetch $t.vcf.gz.tbi data/truth
    # HG002-HG004 high-confidence BEDs exclude sites inconsistent within the trio
    if [ $g = HG005 ]; then fetch $t.bed data/truth; else fetch ${t}_noinconsistent.bed data/truth; fi
  done
  # derived: name-collated FASTQ; secondary and supplementary records dropped,
  # reads without a mate go to a separate file and are not used
  for bam in $d/*_L008.posiSrt.markDup.bam; do
    g=$(basename $bam | sed -E 's/.*Sample_(HG00[0-9]).*/\1/')
    [ -s $d/${g}_R1.fastq.gz ] && continue
    samtools collate -@ 4 -u -O -T $d/tmp_$g $bam \
      | samtools fastq -@ 4 -n -F 0x900 -1 $d/${g}_R1.part.fastq.gz -2 $d/${g}_R2.part.fastq.gz \
          -0 /dev/null -s $d/${g}_singletons.fastq.gz -
    mv $d/${g}_R1.part.fastq.gz $d/${g}_R1.fastq.gz
    mv $d/${g}_R2.part.fastq.gz $d/${g}_R2.fastq.gz
  done
}

extra() {
  # report 03: ClinVar pathogenic sites, difficult-region strata, and the
  # sequencing laboratory's own calls on the HG002-HG005 reads
  local cv=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37
  # dated releases move to archive_2.0/<year>/ after a few weeks
  curl -sSfI $cv/clinvar_20261019.vcf.gz > /dev/null || cv=$cv/archive_2.0/2026
  fetch $cv/clinvar_20261019.vcf.gz data/reference
  fetch $cv/clinvar_20261019.vcf.gz.tbi data/reference
  local st=$giab/release/genome-stratifications/v3.6/GRCh37@all f
  for f in Union/GRCh37_alldifficultregions.bed.gz Union/GRCh37_notinalldifficultregions.bed.gz \
           Union/GRCh37_alllowmapandsegdupregions.bed.gz \
           LowComplexity/GRCh37_AllHomopolymers_ge7bp_imperfectge11bp_slop5.bed.gz \
           SegmentalDuplications/GRCh37_segdups.bed.gz GCcontent/GRCh37_gclt25orgt65_slop50.bed.gz; do
    fetch $st/$f data/reference/stratifications
  done
  local a=$giab_s3/data/AshkenazimTrio/analysis/OsloUniversityHospital_Exome_GATK_jointVC_11242015
  fetch $a/HG002-HG003-HG004.jointVC.filter.vcf data/raw/oslo
  fetch ${a/AshkenazimTrio/ChineseTrio}/HG005.singleSampleVC.filter.vcf data/raw/oslo
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

dbsnp() {
  # GATK b37 resource bundle (contig names as hs37d5). NCBI's own GRCh37 VCFs
  # are 25 GB with RefSeq contig names and came down at under 1 MB/s.
  local b=https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0
  fetch $b/dbsnp_138.b37.vcf.gz data/reference
  fetch $b/dbsnp_138.b37.vcf.gz.tbi data/reference
}

"${1:-panel}"
python3 scripts/check_manifest.py
