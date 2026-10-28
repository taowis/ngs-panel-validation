# Risk register: germline panel pipeline

| | |
|---|---|
| Document | RR-01 |
| Version | 1.1 (1.0 2026-10-28; 1.1 adds R18 after report 05) |
| Effective | 2026-10-28 |
| Author | Kaitao Lai |
| Pipeline version covered | v2 |
| Next review | with any change under SOP-01 section 9, after any entry in `docs/nonconformities.md`, or 12 months |

ISO 15189:2022 asks a laboratory to identify risks to patients from its examination processes, act on them and check that the actions work (clauses 5.6 and 8.5). This register does that for the bioinformatics part of the test only: from FASTQ to filtered calls and the checks in SOP-01 section 7. Sample collection, the wet laboratory, variant classification and reporting are outside it.

Each failure mode is scored as it stands now, with the controls in place, in the style of a failure mode and effects analysis.

- Severity (S): 5 a wrong result that could change a patient's care (a missed or false pathogenic variant); 4 a wrong result likely to be caught before a report; 3 a repeat or a delay; 2 a gap in records or reproducibility; 1 negligible.
- Occurrence (O): 5 expected in most runs; 4 in some runs; 3 occasionally; 2 rarely; 1 not expected.
- Detection (D): 1 always caught by a control before results are used; 2 usually caught; 3 may be caught; 4 unlikely to be caught; 5 no current control can catch it.
- Priority = S x O x D. An action is required when the priority is 40 or more, or when S is 5 and D is 4 or 5. Actions are listed at the end; when an action is closed, the row is scored again.

| ID | Failure mode | Effect | Cause | Controls now | Evidence | S | O | D | Priority |
|---|---|---|---|---|---|---|---|---|---|
| R01 | Sample swap or mislabel | Result reported for the wrong person | Samplesheet error, file mix-up | `09_identity.py` compares calls with the expected genome and checks sex (SOP-01 7.1). For a patient there is no truth set, so only the sex check would apply | Report 03 analysis 7 | 5 | 2 | 4 | 40 |
| R02 | Contaminated sample | False heterozygous calls | Cross-sample DNA or reads | Het/hom ratio and variant counts in run QC ranges (SOP-01 7.2). Low-level contamination is not detected | Report 04, limits | 5 | 2 | 3 | 30 |
| R03 | Panel bases below 20x | Variant missed, reported as negative | Capture gaps, GC content, low yield | Every panel base below 20x listed; ClinVar pathogenic sites below 20x listed; another method needed before a negative result (SOP-01 7.3, 7.4) | Reports 01, 03 | 5 | 4 | 1 | 20 |
| R04 | Hard filter removes true variants | Variant missed | Filter rule not suited to the data | Validation against GIAB truth with intervals (SOP-01 8); change control (SOP-01 9). The SOR rule that did this was removed in change 02 | NC-01; reports 01 to 03 | 5 | 2 | 3 | 30 |
| R05 | Indel error at a long homopolymer | False or missed indel | PCR and sequencing slippage | Listed as a limit in reports 01 and 03; no flag on the calls | Report 01 (566 of 793 indel false positives); report 03 | 5 | 3 | 4 | 60 |
| R06 | Reads from a pseudogene placed on the gene | False calls or missed variants in PMS2 | PMS2CL, segmental duplications | Accuracy measured in difficult regions and in PMS2 (report 03); not stated in SOP-01 as a reporting limit | Report 03 (PMS2 pair in HG002 and HG003) | 5 | 3 | 4 | 60 |
| R07 | Variant type outside scope | CNV, exon deletion, structural or mosaic variant missed | The method cannot detect it | Scope stated in SOP-01 section 2 and the README | SOP-01 | 5 | 3 | 5 | 75 |
| R08 | Wrong or corrupted reference or resource file | Systematic errors | Download fault, wrong version | md5, size and version per file in `data/MANIFEST.tsv`; `check_manifest.py` before a run (SOP-01 6.4) | `data/README.md` | 5 | 1 | 1 | 5 |
| R09 | Software differs from the validated version | Calls change without notice | New machine, environment or image | Locked environment and image; test data run and `check_test_run.py` before first use (SOP-01 6.3); CI on every pull request | `tests/README.md`; CI and local Docker runs 28/28 | 4 | 2 | 1 | 8 |
| R10 | Pipeline changed without validation | Performance changes unknown | Edit to a filter, parameter, tool or region file | Change control: plan first, split data, `05_change_control.py`, report (SOP-01 9); one pull request per change | Report 02 | 5 | 2 | 1 | 10 |
| R11 | Wrong capture BED or parameter file | Wrong reportable region, wrong QC | Parameter file error | Run QC ranges are keyed to the capture file and variant counts fall outside them (SOP-01 7.2) | Report 04 | 4 | 2 | 2 | 16 |
| R12 | Results from an incomplete run used | Missing or partial calls | Interrupted run, memory failure | Run used only when the log has the completion line and the trace has no failed task (SOP-01 6.8); restarts recorded | NC-03 | 4 | 2 | 1 | 8 |
| R13 | A pipeline step overwrites an input file | Wrong comparison figures | Working file named like a staged input | Real input names differ from the working names; test data run | NC-04 | 3 | 1 | 3 | 9 |
| R14 | Performance overstated | Test used beyond what the data support | Few truth variants in the panel (53 SNVs for HG001) | 95% Wilson intervals with every figure; five genomes; pooled intervals marked optimistic (SOP-01 8.5) | Reports 01 to 03 | 4 | 4 | 2 | 32 |
| R15 | Run QC ranges wrong | Good samples stopped or bad samples passed | Ranges set from two to four samples per assay | Ranges marked provisional, to be set again after 20 runs per assay | `config/run_qc.yaml` | 3 | 3 | 3 | 27 |
| R16 | Variant description wrong | Wrong HGVS name on a report | Transcript choice, indel normalisation | None: annotation has not been validated | npaac-mapping, S6.2 | 5 | 3 | 5 | 75 |
| R17 | Records lost or a result cannot be reproduced | Cannot show how a result was produced | Untracked edits, raw data deleted | git with pull requests; run log starts with the commit; data re-creatable from the manifest | SOP-01 11 | 2 | 2 | 2 | 8 |
| R18 | A supplier's update changes output without saying so | Calls, filters or QC change unseen | New release, model or database from an outside supplier | SOP-02: same input through both versions, each run twice, every record compared, not only PASS calls | Report 05: DeepVariant 1.10.0 relabelled 595 records, not in the release notes | 4 | 3 | 2 | 24 |

## Actions

| Risk | Action | Status |
|---|---|---|
| R01 | For patient samples, compare an independent genotype (for example a SNP fingerprint taken at sample receipt) with the calls. Needs a laboratory; recorded as a requirement for any clinical use | open, outside this project |
| R05 | Flag indels next to a homopolymer of 8 bp or more in the output, so they are checked before use (report 03, optional change) | open |
| R06 | Write the PMS2 and segmental duplication limits into SOP-01 section 2 as reporting limits, and require another method (for example long-range PCR) for PMS2 | open |
| R07 | None possible in this pipeline. In a laboratory, every report states the variant types the test does not detect | accepted |
| R16 | Validate HGVS descriptions against a reference set before the pipeline produces them for reports | open, planned |
