# ISO 15189 and this project

ISO 15189:2022 sets the requirements for quality and competence in medical laboratories. Australia adopts it unchanged as AS ISO 15189:2023, and NATA accredits pathology laboratories against it together with the NPAAC standards; NATA's application document (effective February 2026) adds its own interpretations. `docs/npaac-mapping.md` covers the NPAAC requirements for massively parallel sequencing. This document does the same for the ISO 15189 clauses.

This project is not a laboratory and is not accredited, so it cannot conform to ISO 15189. Most of the standard concerns the laboratory as an organisation: impartiality, governance, staff, premises, samples, reports, complaints, audits and management review. The table takes the clauses that bear on an examination's bioinformatics and its quality records, and says where the project meets the requirement, where it does so in part, and where it would need a laboratory. Requirements are summarised in our own words; the standard itself has the text. Clause numbers and titles were checked against NATA's application document and SADCAS's ISO 15189:2022 gap analysis.

## Examination processes (clause 7.3)

| Clause | Requirement, in short | In this project | Status |
|---|---|---|---|
| 7.3.1 | Methods selected for their intended use, with documented criteria where more than one method could be used | Design choices and the alternatives considered: `docs/design-decisions.md`; scope in SOP-01 section 2 | covered |
| 7.3.2 | A validated method verified before use, to show the laboratory obtains the stated performance | Before first use on a new machine, environment or image, the test data are run and compared with expected results (SOP-01 section 6.3); CI repeats this under Docker and Apptainer on every pull request | covered |
| 7.3.3 | A laboratory-developed method validated for its intended use, with performance specifications, and the validation reviewed and recorded | Plans merged before any run; reports 01 to 03 against criteria fixed in `config/acceptance.yaml`. Review by a second person is not possible here: the author filled both roles and the reports say so | partly |
| 7.3.4 | Measurement uncertainty evaluated where relevant, compared with performance specifications and reviewed; where it cannot be evaluated, the reason documented | See "Measurement uncertainty" below | partly |
| 7.3.5 | Reference intervals and clinical decision points defined | Not applicable to germline variant detection. The run QC ranges in `config/run_qc.yaml` are for sample QC, not patient results | not applicable |
| 7.3.6 | Examination procedures documented and available to staff | SOP-01, with document number, version, effective date and review date | covered |
| 7.3.7.1 | Results monitored so that trends and shifts can be seen | Run QC table per sample (report 04). Trend review needs routine runs, which this project does not have | partly |
| 7.3.7.2 | Internal quality control with suitable material, acceptance rules, and results held back when QC fails | Per-sample checks with stop rules: identity, run QC ranges, coverage (SOP-01 section 7). No reference sample is run alongside each batch; a laboratory would include one, for example a GIAB sample at intervals | partly |
| 7.3.7.3 | Participation in external quality assessment, or another comparison where none exists | Not possible outside a laboratory. The comparison with the Oslo laboratory's own calls on the same reads (report 03 analysis 8) is an inter-laboratory comparison but not an EQA programme | not covered |
| 7.3.7.4 | Comparability of results where the same examination is done by different methods or on different equipment | Two libraries of HG001; conda on macOS, Docker and Apptainer on Linux give identical calls on the test data (`tests/README.md`); this pipeline against the Oslo pipeline (report 03) | covered |

## Measurement uncertainty

The test reports a qualitative result (a variant and its genotype at a position), so uncertainty is not a ± value around a number. It is given as follows.

- Every sensitivity, precision and concordance figure has a 95% Wilson interval, and the reports say when a point estimate meets a threshold but its interval does not (SOP-01 section 8.5). With 53 panel SNVs in HG001 even a perfect result has a lower bound of about 0.93.
- Repeatability: the same FASTQ run twice gives identical calls. Reproducibility: two libraries of HG001; four genomes from a second laboratory (reports 01 and 03).
- Dependence on depth: HG002 downsampled to 20x, 30x, 50x and 75x (report 03 analysis 5).
- Dependence on sequence context: accuracy inside and outside GIAB difficult regions, and in PMS2 (report 03 analysis 6).

The 2022 edition also asks for sources outside the analytical step, such as biological variation and pre-examination factors. Here the input is cell-line DNA already sequenced, so DNA quality, extraction and library preparation cannot be assessed. Ancestry is one known source: all samples are of European, Ashkenazi or East Asian ancestry (report 04, limits).

## Data, records and documents

| Clause | Requirement, in short | In this project | Status |
|---|---|---|---|
| 7.6 | Information systems validated before use and after changes; data protected against loss, unauthorised change and access | Every pipeline change validated under SOP-01 section 9; locked software; md5 of every input in `data/MANIFEST.tsv`; git history. Access control, security and backup are not addressed: the data are public cell-line data | partly |
| 7.8 | Plans for continuity of service in an emergency | The pipeline runs on another machine from the image and the manifest, and CI shows it does. There is no written continuity plan | partly |
| 8.3 | Documents identified, reviewed, approved, versioned, and obsolete versions kept from use | SOP-01 and the risk register carry document number, version, effective date and next review; earlier versions are kept in git; every change goes through a pull request. Approval by someone other than the author is not possible here | partly |
| 8.4 | Records created at the time, kept unaltered, retained and retrievable | Reports are not rewritten: report 03 corrects report 01 rather than editing it. Run records with interruptions (`run_record.md`); the first line of every run log is the commit. NATA ties retention times to the NPAAC requirements; raw data here are re-creatable from the manifest rather than kept | covered |

## Risk, nonconforming work and improvement

| Clause | Requirement, in short | In this project | Status |
|---|---|---|---|
| 5.6, 8.5 | Risks to patients identified, acted on in proportion, and the actions checked | `docs/risk-register.md`: 18 failure modes of the bioinformatics part, scored for severity, occurrence and detection, with the controls and the actions still open. Risks in the rest of the laboratory are not covered | partly |
| 7.5 | Nonconforming work identified, its effect assessed, results held back or recalled, and work restarted only once corrected | `docs/nonconformities.md`; SOP-01 section 7 stops a sample on a failed check; section 10 records deviations. No results were ever released, so recall has not arisen | covered |
| 8.6 | Opportunities for improvement identified and acted on | Change 02 (report 02), run QC (report 04), containers and test data | covered |
| 8.7 | Nonconformities investigated for cause, corrected, and the correction's effect reviewed | Four entries in `docs/nonconformities.md`, each with cause, corrective action and how its effect was checked; two closed, two open | covered |

## Externally provided products and services (clause 6.8)

| Clause | Requirement, in short | In this project | Status |
|---|---|---|---|
| 6.8 | Externally provided products and services that affect results are selected, checked and monitored before and during use | SOP-02 verifies a supplier's update from the exported VCFs before use; report 05 applies it to DeepVariant 1.9.0 to 1.10.0 and finds an undocumented change in filter labels. Supplier selection, agreements and performance review need a laboratory | partly |

## Clauses that need a laboratory

Clause 4 (impartiality, confidentiality, requirements regarding patients); clause 5 except 5.6 (legal entity, laboratory director, scope of accreditation, objectives and policies, organisation); clause 6 (personnel and competence, facilities, equipment calibration, reagents, service agreements, and the rest of 6.8), although reference files are handled much like reagent lots, with version and checksum recorded and a check before use; clause 7.2 (pre-examination), 7.4 (reporting and post-examination) and 7.7 (complaints); clause 8.1, 8.2, 8.8 (internal audits and quality indicators) and 8.9 (management review).

## References

International Organization for Standardization. ISO 15189:2022 Medical laboratories - Requirements for quality and competence. 4th ed. Geneva: ISO; 2022. https://www.iso.org/standard/76677.html

Standards Australia. AS ISO 15189:2023 Medical laboratories - Requirements for quality and competence (identical adoption of ISO 15189:2022). Sydney: Standards Australia; 2023. https://store.standards.org.au/product/as-iso-15189-2023 (accessed 2026-10-28).

National Association of Testing Authorities, Australia. General Accreditation Criteria: ISO 15189:2022 Standard Application Document. Issued December 2025, effective February 2026. https://nata.com.au/files/2025/11/ISO-15189-2022-Standard-Application-Document-Effective-Feb-2026.pdf (accessed 2026-10-28).

Southern African Development Community Accreditation Service. SADCAS F 134(d) ISO 15189:2022 Gap Analysis and Transition Plan, Issue 1. 2023-11-18. https://sadcas.org/sites/default/files/2024-03/SADCAS_F_134%28d%29_-_ISO_15189_-_2022_Gap_Analysis_and_Transition_Plan_%5BIssue_1%5D.pdf (accessed 2026-10-28).
