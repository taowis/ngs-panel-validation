# ngs-panel-validation

Before a clinical genetics laboratory reports DNA variants from a sequencing pipeline, it has to show two things: that the pipeline finds the variants that are really there, and that it does not report variants that are not. It has to show this again whenever the pipeline changes, even for a small edit to a filter or a tool version. This project does that work for a small germline pipeline and writes it up the way a diagnostic laboratory would, following the AMP/CAP guidelines for validating NGS bioinformatics pipelines (Roy et al. 2018, reference below).

## What the project does

The pipeline looks for inherited variants in 31 genes linked to hereditary breast, ovarian and colorectal cancer, including BRCA1, BRCA2, PALB2, the mismatch repair genes MLH1, MSH2, MSH6 and PMS2, APC and TP53. The gene list comes from PanelApp Australia.

To score a pipeline you need samples where the right answer is already known. The Genome in a Bottle consortium (GIAB, run by NIST) has sequenced a few reference cell lines many times with different technologies and published a benchmark set of their variants. This project uses HG001 (NA12878), whose exome was sequenced at the Garvan Institute as two separate libraries.

The steps are:

1. Align the sequencing reads to the human reference genome (BWA-MEM).
2. Call variants and apply quality filters (GATK HaplotypeCaller and hard filters).
3. Compare the calls with the GIAB benchmark inside the panel genes (RTG vcfeval). Sensitivity is the share of true variants found; precision is the share of reported variants that are true.
4. Check coverage, and run the same data twice and two libraries of the same person, to see whether results repeat.
5. Judge the results against acceptance criteria written in `config/acceptance.yaml` before any results were produced, so the bar cannot be moved after seeing them.

## What was found

The first version failed. Working out why was the most useful part. All three panel variants missed in the first library had been called correctly and then removed by one rule of the standard GATK filter, a strand bias test (SOR) that misbehaves at homozygous sites. That is a fault in the pipeline. The other failures came from the data rather than the pipeline: the exome was sequenced to about 50x, far shallower than a diagnostic panel, the exome kit does not capture 19% of the panel's coding sequence, and most false indels sit next to long single-base repeats, a known error of PCR-amplified libraries.

The fix went through change control. A written plan came first. Candidate filters were compared on the odd-numbered chromosomes, and the chosen one was tested on the even-numbered chromosomes, so the test data played no part in the choice. Removing the SOR rule raised panel SNV sensitivity from 50/53 and 51/53 (the two libraries) to 53/53 in both, with no loss of precision, and the change was accepted. The remaining failures need data sequenced at diagnostic depth, not pipeline changes. HG001 also has only 53 benchmark variants in these genes, too few for even 53/53 to show the 0.99 target on its own. Report 03 is now testing the pipeline on four more GIAB genomes (HG002 to HG005) sequenced by a second laboratory with a different capture kit.

The rest of this README gives the figures, the layout and how to run it.

## Results in figures

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
- **Documentation**: an SOP (`docs/sop/SOP-01-germline-panel-pipeline.md`) covering running, checks, validation and change control, with templates for plans and reports.

Acceptance criteria are fixed in `config/acceptance.yaml` before any results are generated.

## Layout

```
config/      pipeline parameters, panel definition, acceptance criteria
workflow/    Nextflow pipeline
scripts/     panel BED construction, benchmarking and report scripts
data/        panel BED files (committed); raw data and references (not committed, see data/README.md)
docs/sop/    standard operating procedure and report templates
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

Reports 01 (initial validation) and 02 (SNV filter change) are complete, and SOP-01 is written. Report 03 (HG002 to HG005 from a second laboratory and capture kit) is in progress; its plan is in `docs/reports/03-more-reference-samples/plan.md`.

## Reference

Roy S, et al. Standards and guidelines for validating next-generation sequencing bioinformatics pipelines: a joint recommendation of the Association for Molecular Pathology and the College of American Pathologists. J Mol Diagn. 2018;20(1):4-27. doi:10.1016/j.jmoldx.2017.11.003. PMID: 29154853.
