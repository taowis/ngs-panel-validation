# Report 01: initial validation

Pipeline version: commit 50fc13e (BWA-MEM 0.7.19, GATK 4.6.2.0 HaplotypeCaller, GATK generic hard filters). Run on 2026-10-23 on a 10-core, 16 GB Apple M4, local executor. The run was stopped once and resumed with `-resume`; the run logs are `results/pipeline_info/run_20261023_1026.log` and `run.log`. Criteria: `config/acceptance.yaml`, fixed before the run.

Data: GIAB HG001, Garvan Nextera Rapid Capture Expanded Exome, HiSeq 2500. Two libraries (NIST7035, NIST7086), two lanes each. NIST7035 was processed twice as separate runs. Truth: NIST v4.2.1, GRCh37.

Tables in this folder are produced by `scripts/02_report.py` (accuracy, coverage, reproducibility, acceptance checks) and `scripts/03_diagnose.py` (diagnosis files).

## Outcome

The pipeline does not meet the acceptance criteria. Within-run reproducibility and SNV precision and genotype concordance pass. Sensitivity, indel precision, coverage and between-library concordance fail.

| Criterion | NIST7035 | NIST7086 | Threshold | Result |
|---|---|---|---|---|
| Panel SNV sensitivity | 50/53 = 0.943 (0.846-0.981) | 51/53 = 0.962 (0.873-0.990) | 0.99 | fail |
| Panel SNV precision | 50/50 = 1.000 (0.929-1.000) | 51/51 = 1.000 (0.930-1.000) | 0.99 | pass |
| Panel SNV genotype concordance | 1.000 | 1.000 | 0.99 | pass |
| Panel indels | no truth indels in the panel region | | | not assessable |
| Exome SNV sensitivity / precision | 0.844 / 0.983 | 0.849 / 0.983 | 0.99 | fail |
| Exome SNV genotype concordance | 0.994 | 0.995 | 0.99 | pass |
| Exome indel sensitivity / precision | 0.754 / 0.768 | 0.768 / 0.771 | 0.95 | fail |
| Exome indel genotype concordance | 0.917 | 0.916 | 0.99 | fail |
| Panel bases at >= 20x | 0.845 | 0.861 | 0.98 | fail |
| Within-run, NIST7035 run twice | 102,862 records each, identical md5 | | identical | pass |
| Between-library concordance, exome / panel | 0.896 / 0.944 | | 0.99 | fail |

Intervals are 95% Wilson intervals. NIST7035 run 2 gave the same figures as run 1 throughout.

## Why it fails

The failures have three separate causes. Only the first is a property of the pipeline.

The SNP hard filter removes true calls. All three panel SNVs missed in NIST7035 (5:80168937, 17:41245237 in BRCA1, 17:59763465 in BRIP1) were called with the correct genotype and filtered, each on SOR > 3 alone. Exome-wide, 3,354 of the 7,003 missed SNVs (allele-level) were called and filtered, and 2,796 of those failed only the SOR rule. 2,551 of these were homozygous ALT. At a homozygous site there are no reference reads, so SOR reflects only the strand balance of the ALT reads, which in this capture is often one-sided. See `diagnosis_false_negatives.tsv` and `diagnosis_snv_filter_rules.tsv`.

The data are shallow for this purpose. Mean depth over the capture targets is 47x (NIST7035) and 51x (NIST7086), well below the several hundred-fold depth of a diagnostic panel. 3,426 of the missed SNVs and 811 of the 881 missed indels were not called at all and sit below 20x. The coverage figure is further limited by the capture design: 19% of the panel's padded coding sequence (24,477 of 126,839 bp) is not targeted by this exome kit. MSH2, EPCAM, CHEK2, TP53 and AXIN2 have the lowest fraction at 20x (0.59 to 0.74). Every base below 20x is listed in `low_coverage_regions.tsv`. The between-library comparison is already restricted to positions at 20x or more in both libraries; its discordant calls have not been broken down further here.

Indel false positives are mostly homopolymer errors. Of 793 false-positive indels in NIST7035, 566 are 1 bp insertions or deletions next to a homopolymer of 8 bp or longer (`diagnosis_indel_false_positives.tsv`). This is a known error mode of PCR-amplified short-read libraries. A diagnostic service would normally mask such sites or confirm them by another method.

## Limits of this validation

The panel region holds 53 truth SNVs and no truth indels for HG001. Even 53/53 would give a lower 95% bound of about 0.93, so the panel figures cannot by themselves show sensitivity of 0.99, and panel indel performance cannot be measured. A real validation would add samples (other GIAB genomes, in-house samples with orthogonal results) to reach a few hundred variants per class in the reportable range.

## Follow-up

Change 02 tests removing the SOR rule from the SNP filter, with candidates chosen on odd autosomes and tested on even autosomes (`docs/reports/02-snv-filter-change/`). The depth and capture limits cannot be fixed in the pipeline. They would need panel data sequenced at diagnostic depth.
