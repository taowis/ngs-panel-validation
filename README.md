# NGSPanelValidation

Validation and verification of a germline NGS panel pipeline against the GIAB HG001 (NA12878) reference sample, set out the way a diagnostic laboratory would document it before a pipeline goes live or after a change.

The pipeline itself is ordinary (BWA-MEM, GATK HaplotypeCaller). The point of the project is the evidence around it:

- **Accuracy**: SNV and indel sensitivity and precision against the NIST v4.2.1 benchmark, restricted to the genes of a hereditary cancer panel, with hap.py.
- **Reproducibility**: the two Garvan HG001 exome libraries (NIST7035 and NIST7086) run through the same pipeline, and each library run twice, to show between-run and within-run agreement.
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

## Status

In progress. See `docs/reports/` for completed reports.

## Reference

Roy S, et al. Standards and guidelines for validating next-generation sequencing bioinformatics pipelines: a joint recommendation of the Association for Molecular Pathology and the College of American Pathologists. J Mol Diagn. 2018;20(1):4-27. doi:10.1016/j.jmoldx.2017.11.003. PMID: 29154853.
