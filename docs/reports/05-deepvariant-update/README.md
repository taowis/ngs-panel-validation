# Report 05: verifying a supplier's update, DeepVariant 1.9.0 to 1.10.0

Change: the variant caller supplied by an outside party moves from DeepVariant 1.9.0 to 1.10.0 (WES model). Everything else is held fixed: the same HG003 BAM from pipeline v2, cut to the SureSelect V5 capture; the same reference; the same regions.
Before: `docker://google/deepvariant:1.9.0`, digest `sha256:07e95b34…c2d6ad`. After: `:1.10.0`, digest `sha256:962e5a83…ecbe6`. Full digests in `data/MANIFEST.tsv`.
Procedure: SOP-02, first use. Evaluation: `scripts/12_external_update.py`. Runs: NCI Gadi, 2026-10-28 (`run_record.md`).

Decision: accepted. No call changed. The update relabels 595 non-variant records, which the release notes do not mention; this matters only to software downstream that reads the FILTER column, and is recorded below.

## Why DeepVariant, and why HG003

Laboratories increasingly run variant calling on platforms they did not write (SOPHiA DDM, Illumina DRAGEN and others) and receive new versions from the supplier. Those platforms need a licence, so this report uses a caller that behaves the same way from the laboratory's side: DeepVariant is distributed as versioned container images with a trained model inside, and the laboratory cannot see or change how a call is made. SOP-02 works only from the exported VCFs, so it applies unchanged to a commercial platform.

DeepVariant's WES models for 1.9 and 1.10 were trained on HG001, HG002, HG004, HG005, HG006 and HG007, never on HG003, and never on chromosomes 20 to 22 (training data page in `data/README.md`). Scoring either version on HG002 would therefore flatter it. HG003 is the one genome in this project that neither version has seen.

## Results

`accuracy.tsv`. PASS calls, vcfeval genotype mode; 95% Wilson intervals. The project's own pipeline v2 (GATK HaplotypeCaller with hard filters, report 03) on the same BAM is shown for comparison; it is not part of the decision.

| Region | Type | DeepVariant 1.9.0 | DeepVariant 1.10.0 | GATK v2 |
|---|---|---|---|---|
| Panel | SNV sensitivity | 72/72 = 1.000 (0.949-1.000) | same | 72/72 = 1.000 |
| Panel | SNV precision | 0.986 (1 FP) | same | 0.973 (2 FP) |
| Panel | indel sensitivity / precision | 3/3, 0 FP | same | 2/3, 2 FP |
| Exome | SNV sensitivity / precision | 0.9897 / 0.9970 | same | 0.9820 / 0.9933 |
| Exome | indel sensitivity / precision | 0.9145 / 0.9649 | same | 0.8941 / 0.8945 |
| Exome | SNV / indel genotype concordance | 0.9987 / 0.9753 | same | 0.9985 / 0.9552 |

The two versions give the same counts in every cell. The GATK figures match report 03's HG003 row exactly, which checks that the evaluation here is the one used before.

## Change-control rules

`change_control.tsv`. Every rule passes.

1. No criterion that passed with 1.9.0 fails with 1.10.0.
2. No sensitivity or precision falls by more than 0.005 (none moves at all).
3. Run-to-run noise: each version was run twice, on different Gadi nodes. Each pair of VCFs is identical record for record (md5 of the records in `run_record.md`), so any difference between versions is a real one.
4. PASS calls in the panel targets that differ between versions: none (`panel_discordant_calls.tsv` is empty; 82 PASS records in each).
5. All records compared, PASS or not (`record_differences.tsv`): 101,297 records in each file, the same positions and alleles, and identical QUAL, genotype and every FORMAT value. 595 records differ only in FILTER:

| FILTER in 1.9.0 | FILTER in 1.10.0 | Records | Genotype |
|---|---|---|---|
| RefCall | NoCall | 581 | `./.` |
| PASS | RefCall | 14 | `0/0` |

Rule 5 was added to the script during this report, after the first evaluation showed one fewer PASS record in the exome region with every accuracy figure unchanged. The accuracy figures could not show these changes: vcfeval ignores records with a reference or missing genotype.

## What changed, and does it matter

The 14 records were homozygous reference calls that 1.9.0 marked PASS; 1.10.0 marks them RefCall, which is consistent. None is in the panel targets. A downstream step that takes "FILTER = PASS" as "variant present" without reading the genotype would have reported these 14 sites with 1.9.0. That is a reason to check the genotype as well as the filter, whatever the version.

The 581 records have no genotype (`./.`). 1.10.0 now labels them NoCall. Both versions' VCF headers describe NoCall as "Site has depth=0 resulting in no call", yet these sites have a read depth of 2 to 972 (median 53). A laboratory that counted NoCall sites as coverage failures would now see 581 apparent zero-depth sites per exome. Under SOP-02 this would be raised with the supplier and any QC rule that reads NoCall checked before release.

The release notes for 1.10.0 (https://github.com/google/deepvariant/releases/tag/v1.10.0, read 2026-10-28) describe long-read phasing, ONT, RNA-seq, DeepSomatic and pangenome changes, and say the new multiallelic post-processing is enabled for all modes except WES. They say nothing about FILTER labels.

Why the calls are identical: the WES model's weights are byte-identical in the two images (`gadi/wes_model_md5.txt`: the `variables` files have the same md5; `saved_model.pb`, `fingerprint.pb` and the example-info file differ, so the model was re-exported but not retrained). For whole-exome data, 1.10.0 is the same model with different post-processing labels.

## Deviations from SOP-02

- SOP-02 was written during this work, so the change plan (section 3) was not merged before the runs. The candidates were fixed in advance (the two latest releases, chosen before any run), there was no selection among candidates, and the acceptance rules are the existing ones in `config/acceptance.yaml`, so no choice could be tuned to the result.
- The release notes were read after the first evaluation rather than before the runs.
- Rule 5 was added to `12_external_update.py` after the first evaluation (above).

## Limits

- One genome, one exome kit, one depth. Panel indels: three truth indels, so the indel intervals are wide.
- An update that changes the model would show real call differences; this one does not, so the report shows the procedure working on a quiet update, not on a difficult one. Report 02 is the example of a change that alters calls.
- HGVS descriptions and classifications are not produced here, so SOP-02 section 5 was not exercised.
- Only the WES model was checked. 1.10.0 changes other models, which this says nothing about.

## Files

`accuracy.tsv`, `change_control.tsv`, `repeatability.tsv`, `panel_discordant_calls.tsv`, `record_differences.tsv`, `run_record.md`, and `gadi/` (PBS output, run information, compressed DeepVariant logs, image and model checksums). The scripts that ran on Gadi are in `scripts/gadi/`.
