# Data

Nothing under `data/raw/`, `data/reference/` or `data/truth/` is committed. All of it is public and can be downloaded again from the links below. Access date for all files: 23 October 2026.

Genome build is GRCh37, to match the Garvan exome capture BED.

## Sequencing data

| File | Description | Source |
|---|---|---|
| `NIST7035_TAAGGCGA_L00{1,2}_R{1,2}_001.fastq.gz` | HG001 exome, library 1, two lanes, about 1.8 GB per file | [GIAB Garvan_NA12878_HG001_HiSeq_Exome](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) |
| `NIST7086_CGTACTAG_L00{1,2}_R{1,2}_001.fastq.gz` | HG001 exome, library 2 (separate DNA vial), two lanes | same |
| `nexterarapidcapture_expandedexome_targetedregions.bed.gz` | Illumina Nextera Rapid Capture Expanded Exome targets, hg19 | same |

Sequenced at the Garvan Institute of Medical Research on a HiSeq 2500.

## Truth set

| File | Description | Source |
|---|---|---|
| `HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz` (+ `.tbi`) | NIST v4.2.1 small variant benchmark, autosomes | [GIAB NISTv4.2.1 GRCh37](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh37/) |
| `HG001_GRCh37_1_22_v4.2.1_benchmark.bed` | High-confidence regions | same |

Stratification BEDs (for difficult regions): [GIAB genome-stratifications v3.6](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/).

## Download

```
mkdir -p data/raw data/truth
base=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab
for lib in NIST7035_TAAGGCGA NIST7086_CGTACTAG; do
  for lane in L001 L002; do
    for r in R1 R2; do
      wget -P data/raw $base/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/${lib}_${lane}_${r}_001.fastq.gz
    done
  done
done
wget -P data/raw $base/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz
for f in HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi HG001_GRCh37_1_22_v4.2.1_benchmark.bed; do
  wget -P data/truth $base/release/NA12878_HG001/NISTv4.2.1/GRCh37/$f
done
```

## Committed files

`data/panels/` holds the panel BED files built by `scripts/`, small enough to commit.
