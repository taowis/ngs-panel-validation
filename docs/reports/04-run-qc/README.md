# Report 04: run-level QC metrics and acceptable ranges

The NPAAC requirements for massively parallel sequencing ask a laboratory to define QC metrics for each run, with acceptable ranges and the action taken when a run falls outside them (S6.6, S6.7). The examples given are reads passing filters, the share aligned, the number of SNVs, Ti/Tv and the share in dbSNP. Until now SOP-01 section 7 checked identity, sex, duplicate rate and coverage, with no ranges. This report measures those metrics on the validation runs and sets ranges from them.

Pipeline v2, results from the runs used in reports 02 and 03. Six samples: the two HG001 libraries (Garvan, Nextera; NIST7035_run1 and NIST7086_run1) and HG002 to HG005 (Oslo, SureSelect V5). NIST7035_run2 is the same FASTQ as run1 and its BAM was not kept. Run 2026-10-26 with `scripts/11_run_qc.py`.

## Metrics

| Metric | How |
|---|---|
| `primary_reads` | `samtools flagstat` on the duplicate-marked BAM, primary records |
| `pct_mapped` | primary mapped / primary |
| `pct_properly_paired` | properly paired / paired in sequencing |
| `pct_duplicates` | `PERCENT_DUPLICATION` from MarkDuplicates |
| `mean_capture_depth` | mosdepth, mean over the capture BED |
| `pass_snvs`, `pass_indels` | PASS calls inside the merged capture BED, multiallelic records split |
| `snv_titv` | transitions / transversions of those SNVs |
| `snv_het_hom` | heterozygous / homozygous ALT SNVs |
| `pct_snvs_in_dbsnp`, `pct_indels_in_dbsnp` | share with a dbSNP 138 record at the same position with the same REF and ALT |

dbSNP 138 is the version in the GATK b37 resource bundle; NCBI's current GRCh37 files are 25 GB (see `data/README.md`). The dbSNP match was checked for HG005: all 60,874 matched SNVs have the called allele among the dbSNP record's alleles. SNV and indel counts and Ti/Tv were checked against `bcftools stats` for HG002 and NIST7035_run1 and agree.

## Results

`run_qc.tsv`.

| Sample | Reads (M) | Mapped % | Properly paired % | Duplicates % | Mean depth | PASS SNVs | PASS indels | Ti/Tv | Het/hom | SNVs in dbSNP % | Indels in dbSNP % |
|---|---|---|---|---|---|---|---|---|---|---|---|
| NIST7035_run1 | 81.6 | 99.88 | 99.13 | 11.3 | 47 | 47,637 | 7,387 | 2.41 | 1.61 | 98.5 | 78.2 |
| NIST7086_run1 | 86.9 | 99.87 | 99.07 | 11.3 | 51 | 47,944 | 7,576 | 2.41 | 1.61 | 98.5 | 78.4 |
| HG002 | 150.4 | 99.72 | 98.98 | 4.9 | 140 | 62,785 | 8,809 | 2.42 | 1.57 | 98.5 | 82.7 |
| HG003 | 128.6 | 99.74 | 99.06 | 4.2 | 120 | 62,749 | 8,651 | 2.41 | 1.56 | 98.2 | 82.9 |
| HG004 | 146.1 | 99.75 | 99.04 | 4.9 | 132 | 63,213 | 8,739 | 2.40 | 1.65 | 98.4 | 83.6 |
| HG005 | 122.0 | 99.76 | 99.11 | 5.7 | 114 | 62,156 | 8,426 | 2.44 | 1.34 | 97.9 | 84.2 |

The two assays differ where expected. The Nextera libraries have about twice the duplicate rate at less than half the depth, and the SureSelect data give about 30% more SNVs, from a larger capture region (70.6 Mb with 50 bp added to each bait, against 62.1 Mb) and more depth at its edges. Within an assay the metrics are close, except the het/hom ratio of HG005, which is lower than the others. That is consistent with ancestry: HG005 is Han Chinese, the others European or Ashkenazi.

## Ranges

The ranges are in `config/run_qc.yaml`, one set per assay, since depth, duplicate rate and variant counts depend on the capture kit and the laboratory. They are reference ranges, so they come from these data; this is different from the acceptance criteria, which were fixed before any results. Each limit is the observed value widened by a stated rule (in the file's header), applied on the side that matters: a lower limit for reads, depth, mapping and dbSNP share, an upper limit for duplicates, both sides for variant counts, Ti/Tv and het/hom. All six samples fall within the ranges, as they must by construction. The check was tested by altering values in a copy of the table, and the altered metrics were flagged.

SOP-01 section 7 now requires the check and says what to do when a metric is out of range.

## Limitations

- Two to four samples per assay, and the two Nextera samples are libraries of one person. The ranges are wide on purpose and are provisional; they should be set again from routine runs, for example once there are 20 per assay.
- All samples are cell-line DNA from European, Ashkenazi and East Asian ancestry. Samples of African ancestry carry more variants and more heterozygous sites, and a correct sample could fall outside the variant count and het/hom ranges.
- dbSNP includes variants from earlier studies of some of these cell lines, HG001 in particular, so the dbSNP share may be higher than for a patient sample.
- A sample inside every range can still be wrong. The ranges catch gross failures (a failed library, a wrong capture BED, a broken filter, contamination that raises heterozygosity); the identity and coverage checks are still needed.
