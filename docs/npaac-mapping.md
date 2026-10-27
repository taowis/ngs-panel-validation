# NPAAC requirements and this project

In Australia, medical pathology laboratories are accredited by NATA against AS ISO 15189 together with the standards of the National Pathology Accreditation Advisory Council (NPAAC). For sequencing-based tests the relevant NPAAC document is *Requirements for Human Medical Genome Testing Utilising Massively Parallel Sequencing Technologies* (First Edition 2017, Tier 4). The Australian Commission on Safety and Quality in Health Care, which now publishes the NPAAC standards, lists it as currently in use; all its pathology standards are under review (page accessed 2026-10-25). The PDF it distributes carries a "DRAFT" watermark.

This project is not an accredited laboratory and cannot comply with these requirements: most of them concern a laboratory's governance, staff, wet laboratory, IT and reporting. The table below takes the bioinformatics and quality standards the project can speak to, and says for each where the work is and whether it is covered, partly covered or not covered. Standard numbers are those of the NPAAC document. The ISO 15189 clauses are taken in turn in `docs/iso15189-mapping.md`.

## Bioinformatics (section 6)

| Standard | Requirement, in short | In this project | Status |
|---|---|---|---|
| S6.1 | Extent of validation in proportion to the risk and complexity of the test | Scope limited to germline SNVs and small indels in 31 genes (README, Scope; SOP-01 section 2) | covered |
| S6.2 | Analytical validity of the whole pipeline established | Reports 01 to 03; secondary analysis only, tertiary (annotation, pathogenicity) not done | partly |
| C6.2(ii), (iii)(b) | Performance measured on data not used to optimise the pipeline | Change 02 chosen on odd autosomes, tested on even autosomes; HG002 to HG005 played no part in it and confirmed it (report 02, report 03 analysis 3) | covered |
| C6.2(iii) | TP, FP, FN; sensitivity; precision; reportable range; reproducibility | vcfeval tables per region and variant type, Wilson intervals, within-run and between-library comparisons | covered |
| C6.2(iii) | True negatives, specificity, limit of detection, reference range | vcfeval does not count true negatives, so specificity is not given; limit of detection applies to low allele fraction testing, which is out of scope | not covered |
| S6.3 | Common sources of error identified and mitigated | `03_diagnose.py`: the SOR filter rule (fixed by change 02), homopolymer indels, capture gaps, depth (report 01) | covered |
| S6.4, C6.4(ii) | Reference material used; GIAB NA12878 calls given as the example | GIAB HG001 (NA12878) and HG002 to HG005, NIST v4.2.1 | covered |
| S6.5, C6.5 | Limitations defined: indel size, poor mapping, poor coverage, repeats and homopolymers | Indels up to 50 bp; low-coverage regions listed per gene; accuracy by difficult-region stratum and PMS2 measured in report 03; not yet written into SOP-01 as reporting limits | partly |
| S6.6, S6.7 | Run-level QC metrics with acceptable ranges (reads passing filters, % aligned, number of SNVs, Ti/Tv, % in dbSNP) and actions when not met | SOP-01 section 7: reads, % mapped and properly paired, duplicates, depth, SNV and indel counts, Ti/Tv, het/hom and % in dbSNP checked against per-assay ranges in `config/run_qc.yaml`, with the action when a metric is out of range (report 04). Ranges rest on two to four samples per assay and are provisional | covered |
| S6.8, C6.8(i) | Every change recorded and verified, QC re-established and performance against the previous version documented and explained | Change control in SOP-01 section 9, `05_change_control.py`, report 02 | covered |
| C6.8(ii) | Revision control with reason, date and person for each change | git history, one pull request per change, change plans and reports | covered |
| C6.8(iii) | All components documented: software, scripts, reference sequences, databases | `environment.yml` and `environment.lock.yml`, `data/MANIFEST.tsv` with md5 and version | covered |
| C6.8(iv) | Date of access recorded for public website data | Access dates in `data/MANIFEST.tsv` and `data/README.md` | covered |

## Other sections

| Standard | Requirement, in short | In this project | Status |
|---|---|---|---|
| S3.1 | Scope of testing defined in advance | README, Scope; SOP-01 section 2; criteria in `config/acceptance.yaml` fixed before any results | covered |
| S4.10 | Depth needed for accurate detection established empirically, with false positive and false negative rates | HG002 downsampled to 20x, 30x, 50x and 75x, with sensitivity and precision at each level (report 03 analysis 5) | covered |
| C7.2(iii) | Repeatability and reproducibility determined | Same FASTQ run twice; two libraries of HG001 | covered |
| C7.1(ii) | Cross-platform confirmation (e.g. Sanger) considered | Noted for panel bases below 20x (SOP-01 section 7); not done | not covered |
| S7.3 | External proficiency testing or inter-laboratory exchange | Comparison with the Oslo laboratory's own calls on the same reads (report 03 analysis 8) is similar in spirit but is not an EQA programme | not covered |
| S7.4, S7.5 | QC recorded for each test and deviations explained | Deviations recorded in reports and run records (SOP-01 section 10) | partly |

Not addressed at all, because they need a laboratory: ethics, consent and incidental findings (section 1), personnel (section 2), wet laboratory procedures (section 4, apart from S4.10), IT infrastructure and data security (section 5), and requesting and reporting (section 8).

## Reference

National Pathology Accreditation Advisory Council. Requirements for Human Medical Genome Testing Utilising Massively Parallel Sequencing Technologies (First Edition 2017). Canberra: Commonwealth of Australia; 2017. ISBN 978-1-76007-290-2 (online). https://www.safetyandquality.gov.au/resources/requirements-human-medical-genome-testing-utilising-massively-parallel-sequencing-technologies-first-edition-2017 (accessed 2026-10-25).
