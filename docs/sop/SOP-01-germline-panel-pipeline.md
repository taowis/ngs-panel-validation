# SOP-01: Germline panel pipeline, running, validation and change control

| | |
|---|---|
| Document | SOP-01 |
| Version | 1.0 |
| Effective | 2026-10-24 |
| Author | Kaitao Lai |
| Pipeline version covered | v2 (SNP filter without the SOR rule, report 02) |
| Next review | on any pipeline change, or 12 months |

This procedure was written for this project. It follows the structure a diagnostic laboratory would use, with reference to the AMP/CAP recommendations (Roy et al. 2018), but it is not the procedure of any accredited laboratory.

## 1. Purpose

To describe how the germline small-variant pipeline in this repository is run, how its output is checked before results are used, how it is validated, and how a change to it is verified and recorded.

## 2. Scope

Germline SNVs and small indels (up to about 50 bp) from Illumina short-read exome or panel data, GRCh37 (hs37d5), reported over the 31-gene hereditary cancer panel in `data/panels/` (PanelApp Australia green genes, coding sequence with 20 bp padding).

Out of scope: copy number variants and exon-level deletions, structural variants, mosaic or low allele fraction variants, somatic testing, variants outside the panel targets, variant classification and clinical reporting.

## 3. Responsibilities

The bioinformatician runs the pipeline, performs the checks in section 7, writes the run or validation record and raises deviations. A second person (a scientist or senior bioinformatician) reviews and signs off validation reports and change decisions. In this project both roles are filled by the author; the reports say so.

## 4. Definitions

- Truth set: the GIAB NIST v4.2.1 benchmark for the genome concerned, GRCh37, with its high-confidence BED.
- Evaluable region: the region where accuracy can be measured. Panel evaluable = panel targets & capture & high-confidence. Exome evaluable = capture & high-confidence.
- Sensitivity: truth variants called with the correct genotype / all truth variants in the region (vcfeval genotype mode, PASS calls).
- Precision: correct calls / all PASS calls in the region.
- Genotype concordance: genotype-mode true positives / allele-mode true positives.
- Within-run reproducibility: the same FASTQ run twice with the same pipeline version.
- Between-library reproducibility: two libraries from the same genome compared at positions at >= 20x in both.
- Change: any edit that can alter calls or QC figures: tool version, parameter, filter, reference or resource file, region file, or workflow logic. Edits to documentation and report scripts that do not touch the pipeline are not changes in this sense, but go through the same pull request process.

## 5. Software and data

- Software versions are fixed by `environment.yml` and recorded in full in `environment.lock.yml`: BWA-MEM 0.7.19, samtools/bcftools/htslib 1.24, GATK 4.6.2.0, mosdepth 0.3.14, RTG Tools 3.13, bedtools 2.31.1. Nextflow runs the workflow (`workflow/main.nf`, `nextflow.config`).
- Every downloaded file has a row in `data/MANIFEST.tsv` with its URL, version, access date, md5 and size. `data/README.md` says how each file was obtained and how derived files are built.
- Acceptance criteria are in `config/acceptance.yaml`. They are fixed before a validation is run and changed only with a written reason in a report.

## 6. Running the pipeline

1. Check out the version to be used from `main` and note the commit (`git rev-parse HEAD`). Do not run from a working tree with uncommitted changes.
2. Create the environment once: `conda env create -f environment.yml`. The environment path must not contain spaces.
3. Fetch or check the inputs: `scripts/00_download.sh <step>`, then `python3 scripts/check_manifest.py`. Any MISMATCH stops the procedure until it is explained.
4. Prepare the samplesheet (`sample,library,lane,fastq_1,fastq_2`, one row per lane). For a genome other than HG001, write a parameter file in `config/params/` naming the samplesheet, capture BED, evaluable region, truth set and output folder (see `config/params/HG002.yaml`).
5. Start the run so that it survives the terminal closing and the commit is the first line of the log:

   ```
   nohup sh -c 'git rev-parse HEAD; caffeinate -i nextflow run workflow/main.nf -profile local -resume \
       -ansi-log false [-params-file config/params/<genome>.yaml]' \
       > <outdir>/pipeline_info/run.log 2>&1 < /dev/null &
   ```

   `caffeinate` applies to macOS; on a server, run under the scheduler instead.
6. If a run is interrupted, restart it with the same command; `-resume` reuses completed tasks. Record the interruption and restart times in the run record.
7. The run is complete when the log ends with the Nextflow completion line and `trace.tsv` shows no failed tasks. A run with a failed task is not used; the cause is recorded and the run repeated.

## 7. Checks before results are used

For every sample:

1. Sample identity and sex: `scripts/09_identity.py` (calls against the expected genome's truth set, chrX heterozygosity, chrY depth). A mismatch stops the procedure; nothing from that sample is used until it is resolved.
2. Duplicate rate from `<outdir>/bam/<sample>.markdup_metrics.txt`. Record it; a rate far outside earlier runs on the same assay is investigated before going on.
3. Coverage: `scripts/02_report.py` writes `coverage_by_gene.tsv` and `low_coverage_regions.tsv`. Every panel base below 20x is listed. On a patient sample these bases would need another method (for example Sanger sequencing) before a negative result could be reported for them.
4. Known pathogenic sites below 20x: `scripts/07_clinvar_sites.py` lists ClinVar pathogenic and likely pathogenic positions in the panel that are not covered.

## 8. Validation

A validation establishes the performance of a pipeline version before it is used, or for a new assay or sample type.

1. Write down the purpose, data, regions, criteria and analyses before running anything. Use `templates/validation-plan.md`. Merge the plan into `main` before the first run.
2. Run the pipeline (section 6) on the reference samples. Within-run reproducibility needs one sample run twice; between-library reproducibility needs two libraries of one genome.
3. Produce the tables with `scripts/02_report.py --name <report> --results <outdir> [...]`. Where criteria fail, run `scripts/03_diagnose.py` to separate pipeline causes (filtering, calling) from data causes (depth, capture, sequence context).
4. Write the report from `templates/validation-report.md`: version and commit, data, deviations from the plan, results against each criterion with 95% Wilson intervals, diagnosis of failures, limits, and the outcome. A failed validation is reported as failed.
5. Give intervals as well as point estimates. With few truth variants in the reportable range (53 SNVs for HG001), a perfect result still has a wide interval; say so rather than claim the threshold is met.

## 9. Change control

1. Raise the change as a plan (`templates/change-plan.md`): the reason, the exact candidates, the data, the selection rule, and the acceptance rules from `config/acceptance.yaml`. Merge the plan before any candidate is run.
2. If candidates are being chosen on the validation data, split the data so that the choice and the test use different variants (report 02 used odd and even autosomes), with `scripts/04_split_eval.py`.
3. Make the selected change on a feature branch, with a comment in the code naming the change record, and merge it through a pull request.
4. Rerun the pipeline on the same data into a new output folder. Keep the earlier output.
5. Run `scripts/05_change_control.py`. Every rule must pass: no criterion that passed before now fails, no sensitivity or precision drop above `max_drop` on the test set, within-run reproducibility unchanged, and every call that differs inside the panel listed.
6. Review each differing panel call against the truth set and the read evidence (depth, allele counts, mapping quality). Record the review in the report.
7. Write the change report (`templates/change-report.md`) with the decision: accepted or rejected. A rejected change is reverted and the report kept.
8. Where later data that played no part in the choice becomes available, check the change on it as well (report 03, analysis 3).

## 10. Deviations

Anything that departs from this SOP or from a plan is recorded in the report concerned, with what happened, why, and what was done: for example an interrupted run, a changed download source, a script fault found after results were produced and the correction made. Results are never quietly regenerated.

## 11. Records

| Record | Where |
|---|---|
| Code, configuration, criteria | git history, `main` branch; every change by pull request |
| Plans and reports | `docs/reports/<NN-name>/` |
| Run logs, trace and Nextflow report | `<outdir>/pipeline_info/` (not committed; summary tables are copied into the report folder) |
| Data sources | `data/MANIFEST.tsv`, `data/README.md` |
| Raw data, BAMs, VCFs | not committed; re-creatable from the manifest and the commit |

## 12. References

Roy S, et al. Standards and guidelines for validating next-generation sequencing bioinformatics pipelines: a joint recommendation of the Association for Molecular Pathology and the College of American Pathologists. J Mol Diagn. 2018;20(1):4-27. doi:10.1016/j.jmoldx.2017.11.003. PMID: 29154853.
