# Report 03: more reference samples, plan

Written 2026-10-24, before any of the samples below were run through the pipeline.

## Reason

Report 02 left two questions open.

1. The panel region holds only 53 truth SNVs and no truth indels for HG001, so even 53/53 gives a lower 95% bound of about 0.93. More genomes are needed to say anything closer to 0.99.
2. Report 01 put the exome-wide failures (SNV and indel sensitivity, coverage) down to the depth and capture of the Garvan data rather than to the pipeline. That was an explanation, not a test. Running the same pipeline on data from another laboratory and capture kit tests it.

A third use comes for free. Change 02 (removing the SOR rule) was chosen and tested on HG001 only. These samples played no part in it, so they give an independent check of that change.

## Samples

GIAB exomes sequenced at Oslo University Hospital in 2015 (`data/README.md`): Agilent SureSelect Human All Exon V5, HiSeq 2500, one lane each, 2 x 126 bp (the GIAB README says 150 bp; the reads in the BAMs are 126 bp).

| Genome | Description | Truth set |
|---|---|---|
| HG002 | Ashkenazi son | NIST v4.2.1 GRCh37, `_noinconsistent` high-confidence BED |
| HG003 | Ashkenazi father | as above |
| HG004 | Ashkenazi mother | as above |
| HG005 | Han Chinese son | NIST v4.2.1 GRCh37 |

GIAB provides only the laboratory's aligned BAMs. The reads are taken back to FASTQ (`scripts/00_download.sh oslo`) and run through the whole pipeline from alignment, as for HG001.

HG001 is represented by NIST7035_run1 from `results_v2/` (commit 26941b8). It is not rerun.

## Pipeline

v2, as accepted in report 02. The filters and tools are unchanged. The only code changes are that the sample name in the read group comes from a parameter instead of being fixed to HG001, and that each genome gets a parameter file naming its samplesheet, capture BED, truth set and evaluable region. With the default parameters the commands for HG001 are the same as in 26941b8.

Each genome is run separately, into `results_v2_<genome>/`.

## Regions

- Capture: the SureSelect V5 bait intervals with 50 bp added each side, as published by GIAB with the Oslo data (the plain target BED is not distributed there). This is wider than the Nextera "targeted regions" BED used for HG001, so the exome region includes more target flank, where depth falls off. HaplotypeCaller adds its usual 100 bp padding.
- Exome evaluable: capture & that genome's high-confidence BED.
- Panel evaluable: panel targets (`data/panels/panel_targets.bed`, unchanged) & capture & high-confidence, built per genome by `01_build_panel.py`. Truth variant counts per gene go in `data/panels/panel_summary.<genome>.tsv`.

## Analyses

Fixed here, before any results:

1. Per genome, every accuracy and coverage criterion in `config/acceptance.yaml`, with the same thresholds. The reproducibility criteria are not assessable, as each genome has one library and one run.
2. Pooled panel accuracy over HG001 (NIST7035_run1) and HG002 to HG005: true positives, false negatives and false positives summed over the five genomes, for SNVs and for indels, with 95% Wilson intervals. The point estimate is judged against the same threshold, as in `acceptance.yaml`.
3. Independent check of change 02: the v1 SNP filter (with the SOR rule) applied to the same calls of HG002 to HG005, compared with v2 over the exome evaluable region and the panel evaluable region.

If any of HG002 to HG005 misses criteria, `03_diagnose.py` is run on it, as for report 01.

## How the results will be read

- Analysis 1: if the exome-wide SNV and indel criteria and the coverage criterion pass on the Oslo data, that supports report 01's explanation (the data, not the pipeline). If they fail, the diagnosis decides whether the cause is again depth and capture, or something in the pipeline that report 01 missed.
- Analysis 2: HG002 inherits one allele at each site from HG003 and one from HG004, so many panel positions are counted more than once in the pool. The pooled figures will be given with the number of distinct positions, and the Wilson interval, which assumes independent variants, will be described as optimistic.
- Analysis 3: change 02 is supported if v2's SNV sensitivity is at least v1's in every genome and its SNV precision is no more than 0.005 below (the `max_drop` in `acceptance.yaml`). If precision falls by more than that in any genome, the result is reported as it is and change 02 is reopened as a new change; v2 is not edited quietly.

The outcome is written up in this folder whatever it is.
