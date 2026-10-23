# Report 02: SNV hard filter change

Change: remove the SOR > 3 rule from the SNP hard filter (indel filter unchanged).
Before: v1, commit 50fc13e, results in `results/`, report 01.
After: v2, commit 26941b8, results in `results_v2/`, run 2026-10-23 on the same data and machine. Alignment, duplicate marking, variant calling and coverage were reused from the v1 run (same inputs and code, Nextflow cache). Filtering and everything after it were rerun.

The plan (`plan.md`) was merged before any candidate was run (PR #6). Decision: **accepted**.

## Choosing the candidate (development set)

Odd autosomes, exome evaluable region minus panel targets, NIST7035_run1 (`dev_candidates.tsv`):

| Candidate | SNV sensitivity | SNV precision | F |
|---|---|---|---|
| v1 | 0.8404 | 0.9820 | 0.9057 |
| A, SOR not applied to homozygous ALT | 0.8959 | 0.9826 | 0.9373 |
| B, no SOR rule | 0.9005 | 0.9809 | 0.9389 |

B has the higher F-measure and its precision is 0.0011 below v1, inside the 0.005 allowed, so B was selected as the plan required. Indel results were the same for all three candidates.

## Test set

Even autosomes, same region definition, all runs, pipeline outputs as produced (`test_pipeline.tsv`):

| Run | Type | Sensitivity v1 → v2 | Precision v1 → v2 |
|---|---|---|---|
| NIST7035_run1 | SNV | 0.8471 → 0.9064 | 0.9833 → 0.9819 |
| NIST7086_run1 | SNV | 0.8515 → 0.9115 | 0.9838 → 0.9821 |
| NIST7035_run1 | indel | 0.7571 → 0.7571 | 0.7666 → 0.7666 |
| NIST7086_run1 | indel | 0.7738 → 0.7738 | 0.7705 → 0.7705 |

NIST7035_run2 is identical to run 1. On held-out chromosomes, SNV sensitivity rises by about 6 points. SNV precision falls by 0.0014 to 0.0017: about 50 more false positives against about 1,250 more true positives per run.

## Change-control rules

All rules in `config/acceptance.yaml` hold (`change_control.tsv`):

1. No criterion that passed in v1 fails in v2. Ten checks, all still pass.
2. No sensitivity or precision drop above 0.005 on the test set. The largest drop is 0.0017 (SNV precision, NIST7086).
3. Every PASS call that differs inside the panel targets is listed in `panel_discordant_calls.tsv` and reviewed below.
4. Within-run reproducibility: the two NIST7035 runs still give identical VCF records.

## Review of the panel calls that changed

There are 19 rows covering 7 positions. Every change goes from filtered (SNP_HARD) in v1 to PASS in v2. No call was lost and no genotype changed.

| Position | Gene | Genotype called | Truth | Note |
|---|---|---|---|---|
| 5:80168937 G>A | MSH3 | 1/1 | 1/1 | missed in v1, panel evaluable region |
| 17:41245237 A>G | BRCA1 | 0/1 | 0/1 | missed in v1 (NIST7035), panel evaluable region |
| 17:59763465 T>C | BRIP1 | 1/1 | 1/1 | missed in v1, panel evaluable region |
| 12:133208886 C>T | POLE | 1/1 | 1/1 | panel target, not in the capture design |
| 12:133233705 T>C | POLE | 1/1 | 1/1 | panel target, not in the capture design |
| 13:32936646 T>C | BRCA2 | 1/1 | 1/1 | panel target, not in the capture design |
| 17:29553485 G>A | NF1 | 1/1 | outside high-confidence | NIST7035 only. NIST7086 calls it 1/1 PASS in both versions (43 ALT reads, no REF reads) |

Six of the seven agree with the truth genotype. The NF1 call cannot be scored against the benchmark. It is a confident homozygous call (QD 31, MQ 60, 40 ALT reads and no REF reads), and the other library makes the same call, so there is no evidence that it is an error. On a live test it would go to the usual review of calls outside validated regions.

## Results of v2 against the acceptance criteria

The full tables are in this folder (`accuracy.tsv`, `acceptance_checks.tsv` and others, produced by `02_report.py --results results_v2`).

| Criterion | v1 | v2 | Threshold |
|---|---|---|---|
| Panel SNV sensitivity, NIST7035 / NIST7086 | 0.943 / 0.962 | 53/53 = 1.000 (0.932-1.000), both | 0.99 |
| Panel SNV precision | 1.000 | 1.000 | 0.99 |
| Exome SNV sensitivity, NIST7035 | 0.844 | 0.903 | 0.99 |
| Exome SNV precision, NIST7035 | 0.983 | 0.981 | 0.99 |
| Exome indel sensitivity / precision | 0.754 / 0.768 | unchanged | 0.95 |
| Panel bases at >= 20x | 0.845 / 0.861 | unchanged | 0.98 |
| Between-library concordance, exome / panel | 0.896 / 0.944 | 0.927 / 0.982 | 0.99 |

The panel SNV criteria now pass on the point estimate. With only 53 truth SNVs, the lower confidence bound (0.93) still does not show sensitivity of 0.99. Exome-wide sensitivity, indel performance, coverage and between-library concordance still fail. Report 01 explains these as the depth and capture of this dataset and homopolymer indel errors. The change does not address them, so v2 is still not validated for clinical use on this data.

## Decision

The change is accepted as an improvement over v1 under the change-control rules, and v2 becomes the pipeline default. This is not a claim that v2 meets the acceptance criteria.
