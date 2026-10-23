# ngs-panel-validation

Validation and verification of a germline NGS panel pipeline against the GIAB HG001 (NA12878) reference sample, set out the way a diagnostic laboratory would document it before a pipeline goes live or after a change.

## Results so far

GIAB HG001, Garvan Nextera exome (about 50x), 31-gene hereditary cancer panel from PanelApp Australia. Figures are for NIST7035 where only one is given; full figures and reasoning are in `docs/reports/`.

| | v1 (report 01) | v2 (report 02) | Criterion |
|---|---|---|---|
| Panel SNV sensitivity | 50/53, 51/53 | 53/53 (95% CI 0.93-1.00) | 0.99 |
| Panel SNV precision | 1.000 | 1.000 | 0.99 |
| Exome SNV sensitivity / precision | 0.844 / 0.983 | 0.903 / 0.981 | 0.99 |
| Exome indel sensitivity / precision | 0.754 / 0.768 | unchanged | 0.95 |
| Panel bases at >= 20x | 0.845 | unchanged | 0.98 |
| Within-run (same data run twice) | identical | identical | identical |
| Between-library concordance, panel | 0.944 | 0.982 | 0.99 |

v1 failed validation. The diagnosis separated a pipeline fault (the SOR rule of the generic SNP hard filter removed true homozygous calls) from limits of the data (shallow exome, 19% of the panel's coding sequence not captured, homopolymer indel errors). Change 02 removed the SOR rule. The candidate was chosen on odd autosomes and tested on even autosomes against a plan fixed beforehand, and it passed change control. The remaining failures need data at diagnostic depth, not pipeline changes.

## Approach

The pipeline itself is ordinary (BWA-MEM, GATK HaplotypeCaller). The point of the project is the evidence around it:

- **Accuracy**: SNV and indel sensitivity and precision against the NIST v4.2.1 benchmark, restricted to the genes of a hereditary cancer panel, with RTG vcfeval.
- **Reproducibility**: the two Garvan HG001 exome libraries (NIST7035 and NIST7086) compared with each other, and NIST7035 run twice end to end, to show between-library and within-run agreement.
- **Change control**: a pipeline change (tool version or parameter) run against the same data, with a before/after verification report and an accept or reject decision against predefined criteria.
- **Documentation**: an SOP and report templates in `docs/`.

Acceptance criteria are fixed in `config/acceptance.yaml` before any results are generated.

## Layout

```
config/      pipeline parameters, panel definition, acceptance criteria
workflow/    Nextflow pipeline
scripts/     panel BED construction, benchmarking and report scripts
data/        panel BED files (committed); raw data and references (not committed, see data/README.md)
docs/sop/    standard operating procedure
docs/reports/ validation and verification reports
```

## Running

```
conda env create -f environment.yml          # env path must not contain spaces
scripts/00_download.sh panel && scripts/00_download.sh reference && scripts/00_download.sh reads
python3 scripts/01_build_panel.py
nextflow run workflow/main.nf -profile local
python scripts/02_report.py --name <report-name>
python scripts/03_diagnose.py --name <report-name>    # why calls were missed or false
```

For a pipeline change: write the plan, choose between candidates with `scripts/04_split_eval.py --set dev`, rerun the pipeline with `--outdir results_v2`, then `04_split_eval.py --set test --dirs results results_v2` and `scripts/05_change_control.py`. Run the Python scripts inside the conda env.

`config/samplesheet.csv` lists the lanes; `nextflow.config` holds paths and resource limits.

## Status

Reports 01 (initial validation) and 02 (SNV filter change) are complete. An SOP is still to be written in `docs/sop/`.

## Reference

Roy S, et al. Standards and guidelines for validating next-generation sequencing bioinformatics pipelines: a joint recommendation of the Association for Molecular Pathology and the College of American Pathologists. J Mol Diagn. 2018;20(1):4-27. doi:10.1016/j.jmoldx.2017.11.003. PMID: 29154853.
