# Data

Every downloaded file has a row in `data/MANIFEST.tsv` with its exact URL, version, access date, md5 and size. To fetch everything again and check it:

```
scripts/00_download.sh panel     # annotation, capture BED, truth set (~200 MB)
scripts/00_download.sh reads     # HG001 exome FASTQ (~30 GB, run on HPC)
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
| `NIST7035_TAAGGCGA_L00{1,2}_R{1,2}_001.fastq.gz` | library 1, two lanes, about 1.8 GB per file | not yet downloaded |
| `NIST7086_CGTACTAG_L00{1,2}_R{1,2}_001.fastq.gz` | library 2 (separate DNA vial), two lanes | not yet downloaded |

The FASTQ are fetched from GIAB's AWS mirror, https://giab.s3.amazonaws.com/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/, which serves the same files (identical sizes, checked 2026-10-23) at about 7 MB/s against 0.1-0.2 MB/s from ftp-trace. FASTQ rows go into the manifest when they are downloaded.

## Reference genome: hs37d5 (GRCh37 + decoy)

`data/reference/hs37d5.fa.gz` from the 1000 Genomes AWS mirror, https://1000genomes.s3.amazonaws.com/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz (original 2011 upload, 901,788,727 bytes). EBI's copy at https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz is a 2017 recompression (892,326,179 bytes); both serve the same `hs37d5.fa.gz.fai` (md5 bb77e60e9a492fd0172e2b11e6c16afd), so the uncompressed sequence layout is the same.

`scripts/00_download.sh reference` also builds the derived files: `hs37d5.fa`, `.fai`, `.dict` (GATK), BWA index and `hs37d5.sdf` (RTG).

## Truth set: GIAB HG001 NIST v4.2.1, GRCh37

From https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh37/ (accessed 2026-10-23):

- `HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz` and `.tbi`: small variant benchmark, autosomes
- `HG001_GRCh37_1_22_v4.2.1_benchmark.bed`: high-confidence regions

Stratification BEDs for difficult regions, if needed later: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/

## Committed files built here

`data/panels/panel_genes.tsv`, `panel_targets.bed`, `panel_evaluable.bed` and `panel_summary.tsv` are built by `scripts/01_build_panel.py` from the files above.
