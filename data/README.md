# Data

Every downloaded file has a row in `data/MANIFEST.tsv` with its exact URL, version, access date, md5 and size. To fetch everything again and check it:

```
scripts/00_download.sh panel       # annotation, capture BED, truth set (~200 MB)
scripts/00_download.sh reference   # hs37d5 FASTA (~900 MB), then builds the indexes
scripts/00_download.sh reads       # HG001 exome FASTQ (~16 GB)
scripts/00_download.sh oslo        # HG002-HG005 exome BAMs (~38 GB), FASTQ and truth sets
scripts/00_download.sh extra       # ClinVar, GIAB stratifications, Oslo's own calls
scripts/00_download.sh dbsnp       # dbSNP 138, b37 (~1.5 GB)
python3 scripts/check_manifest.py
```

Nothing under `data/raw/`, `data/reference/` or `data/truth/` is committed. Genome build is GRCh37, to match the Garvan exome capture BED.

## Gene panel: PanelApp Australia

Green (diagnostic grade) genes from three panels, fetched by pinned version on 2026-10-23:

| Panel | ID | Version | URL |
|---|---|---|---|
| Breast Cancer | 4375 | 2.0 | https://panelapp-aus.org/api/v1/panels/4375/?version=2.0 |
| Ovarian Cancer | 4374 | 2.0 | https://panelapp-aus.org/api/v1/panels/4374/?version=2.0 |
| Colorectal Cancer and Polyposis | 4371 | 2.1 | https://panelapp-aus.org/api/v1/panels/4371/?version=2.1 |

The API responses are committed in `data/panels/panelapp/` and are the record; `scripts/00_download.sh panelapp` refetches them.

## Gene annotation: Ensembl GRCh37 release 87

`data/reference/Homo_sapiens.GRCh37.87.gtf.gz` from https://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz (accessed 2026-10-23). Release 87 matches the GRCh37 Ensembl gene IDs that PanelApp reports.

## Sequencing data: GIAB HG001 exome (Garvan)

From https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/. Sequenced at the Garvan Institute of Medical Research on a HiSeq 2500, Nextera Rapid Capture Expanded Exome.

| File | Description | Accessed |
|---|---|---|
| `nexterarapidcapture_expandedexome_targetedregions.bed.gz` | capture targets, hg19 (`chr` prefix stripped when read) | 2026-10-23 |
| `NIST7035_TAAGGCGA_L00{1,2}_R{1,2}_001.fastq.gz` | library 1, two lanes, about 1.8 GB per file | 2026-10-23 |
| `NIST7086_CGTACTAG_L00{1,2}_R{1,2}_001.fastq.gz` | library 2 (separate DNA vial), two lanes | 2026-10-23 |

The FASTQ are fetched from GIAB's AWS mirror, https://giab.s3.amazonaws.com/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/, which serves the same files (identical sizes, checked 2026-10-23) at about 7 MB/s against 0.1-0.2 MB/s from ftp-trace. The md5 and size of each FASTQ are in the manifest.

## Reference genome: hs37d5 (GRCh37 + decoy)

`data/reference/hs37d5.fa.gz` from the 1000 Genomes AWS mirror, https://1000genomes.s3.amazonaws.com/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz (original 2011 upload, 901,788,727 bytes). EBI's copy at https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz is a 2017 recompression (892,326,179 bytes); both serve the same `hs37d5.fa.gz.fai` (md5 bb77e60e9a492fd0172e2b11e6c16afd), so the uncompressed sequence layout is the same.

`scripts/00_download.sh reference` also builds the derived files: `hs37d5.fa`, `.fai`, `.dict` (GATK), BWA index and `hs37d5.sdf` (RTG).

## Truth set: GIAB HG001 NIST v4.2.1, GRCh37

From https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh37/ (accessed 2026-10-23):

- `HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz` and `.tbi`: small variant benchmark, autosomes
- `HG001_GRCh37_1_22_v4.2.1_benchmark.bed`: high-confidence regions

Stratification BEDs for difficult regions, if needed later: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/

## Sequencing data: GIAB HG002 to HG005 exomes (Oslo University Hospital)

Used in report 03. Agilent SureSelect Human All Exon V5, HiSeq 2500, one lane per sample, 2 x 126 bp, sequenced in 2015. GIAB publishes only the laboratory's duplicate-marked BAMs (BWA-MEM 0.7.12 against b37 + decoy), from https://giab.s3.amazonaws.com/data/AshkenazimTrio/<sample>/OsloUniversityHospital_Exome/ and https://giab.s3.amazonaws.com/data/ChineseTrio/HG005_NA24631_son/OsloUniversityHospital_Exome/ (the same paths exist on ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/). `scripts/00_download.sh oslo` downloads them to `data/raw/oslo/`, accessed 2026-10-24, and converts each back to paired FASTQ (`samtools collate | samtools fastq`, secondary and supplementary records dropped). The pipeline realigns the FASTQ. On 2026-10-24 the HG002 and HG005 downloads from S3 slowed to about 0.5 MB/s each; they were stopped and finished as parallel byte-range requests appended to the partial files. All four BAMs match their S3 multipart ETags, and their md5 are in the manifest.

| File | Genome | Size |
|---|---|---|
| `151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam` | HG002, Ashkenazi son | 10.3 GB |
| `..._HG003-EEogPU_v02-KIT-Av5_TCTTCACA_L008.posiSrt.markDup.bam` | HG003, Ashkenazi father | 9.0 GB |
| `..._HG004-EEogPU_v02-KIT-Av5_CCGAAGTA_L008.posiSrt.markDup.bam` | HG004, Ashkenazi mother | 10.2 GB |
| `..._HG005-EEogPU_v02-KIT-Av5_CGCATACA_L008.posiSrt.markDup.bam` | HG005, Han Chinese son | 8.9 GB |

Capture regions: `wex_Agilent_SureSelect_v05_b37.baits.slop50.merged.list` from https://giab.s3.amazonaws.com/data/AshkenazimTrio/analysis/OsloUniversityHospital_Exome_GATK_jointVC_11242015/, a Picard interval list of the V5 baits with 50 bp added each side, as the laboratory used for calling. The plain target BED is not distributed with the data. The script converts it to `sureselect_v5_b37.baits.slop50.bed.gz`.

The laboratory's own calls on the same reads, from the same folder (and the ChineseTrio equivalent): `HG002-HG003-HG004.jointVC.filter.vcf` (joint HaplotypeCaller calls, VQSR) and `HG005.singleSampleVC.filter.vcf`. Used only for the comparison in report 03.

## Truth sets: GIAB HG002 to HG005 NIST v4.2.1, GRCh37

From https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/<sample>/NISTv4.2.1/GRCh37/ and .../ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh37/, via the GIAB S3 mirror (same paths under https://giab.s3.amazonaws.com/release/), accessed 2026-10-24: `<genome>_GRCh37_1_22_v4.2.1_benchmark.vcf.gz` and `.tbi`, and the high-confidence BED. For HG002, HG003 and HG004 this is `_benchmark_noinconsistent.bed`, which leaves out sites where the three benchmarks are inconsistent with Mendelian inheritance; HG005 has `_benchmark.bed`.

## ClinVar, GRCh37

`data/reference/clinvar_20261019.vcf.gz` and `.tbi` from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20261019.vcf.gz (accessed 2026-10-24). ClinVar moves dated releases to https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2026/ after a few weeks. Used for the pathogenic site coverage in report 03.

## dbSNP build 138, b37

`data/reference/dbsnp_138.b37.vcf.gz` and `.tbi` from the GATK resource bundle that the Broad Institute keeps on Google Cloud, https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz (accessed 2026-10-26). Its contig names match hs37d5. Used only for the known-variant fraction in the run QC (`scripts/11_run_qc.py`, report 04). NCBI's own GRCh37 VCFs (https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/GCF_000001405.25.gz) are 25 GB, use RefSeq contig names and downloaded at under 1 MB/s from here, so they were not used.

## GIAB genome stratifications v3.6, GRCh37

From https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh37@all/ (accessed 2026-10-24), to `data/reference/stratifications/`: `Union/GRCh37_alldifficultregions.bed.gz`, `Union/GRCh37_notinalldifficultregions.bed.gz`, `Union/GRCh37_alllowmapandsegdupregions.bed.gz`, `LowComplexity/GRCh37_AllHomopolymers_ge7bp_imperfectge11bp_slop5.bed.gz`, `SegmentalDuplications/GRCh37_segdups.bed.gz`, `GCcontent/GRCh37_gclt25orgt65_slop50.bed.gz`, checked against `GRCh37-genome-stratifications-md5s.txt` from the same folder.

## Committed files built here

`data/panels/panel_genes.tsv`, `panel_targets.bed`, `panel_evaluable.bed` and `panel_summary.tsv` are built by `scripts/01_build_panel.py` from the files above. `panel_evaluable.<genome>.bed` and `panel_summary.<genome>.tsv` are the same for HG002 to HG005, with the SureSelect V5 capture and that genome's truth set (`--tag`; commands in the script's docstring).
