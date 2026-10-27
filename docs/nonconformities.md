# Nonconformities and corrective actions

A nonconformity is anything found to be wrong with the pipeline, its results, its documents or the way a procedure was carried out: a fault in a script, a result that does not meet a requirement, a conclusion later shown to be incomplete, or a run that departed from SOP-01. Deviations that are planned and explained in a report are recorded there and not repeated here. Entries follow `docs/sop/templates/nonconformity.md` and are added, not edited: a later finding goes in a new entry that refers to the earlier one. Each entry is checked against `docs/risk-register.md`.

## NC-01: SNP hard filter removed true variants

| | |
|---|---|
| Found | 2026-10-23, initial validation (report 01, PR #8) |
| What | All three panel SNVs missed in NIST7035 were called with the right genotype and then filtered, each on SOR > 3 alone. Exome-wide the rule caused 3,354 of the 7,003 missed SNVs. |
| Results affected | Pipeline v1 only, which was never released; report 01 failed |
| Immediate action | Pipeline not released; report 01 written as a failed validation |
| Cause | The SOR threshold in the GATK hard filter recommendation removes true calls in this data |
| Corrective action | Change 02 under SOP-01 section 9: plan merged first (PR #6), candidates chosen on odd autosomes and tested on even autosomes, change accepted (report 02, PR #9) |
| Effectiveness | Checked on four genomes that played no part in the change: SNV sensitivity up 4 to 5 points in every genome, precision down at most 0.0003, four panel SNVs recovered (report 03 analysis 3, PR #25) |
| Risk register | R04 |
| Status | closed 2026-10-26 |

## NC-02: Report 01 gave an incomplete cause for the exome-wide failures

| | |
|---|---|
| Found | 2026-10-26, report 03 |
| What | Report 01 put HG001's exome-wide failures down to shallow data, capture gaps and homopolymer indels, and said the remaining failures needed data at diagnostic depth. The Oslo data at 114x to 140x still fail, and the failures sit in segmental duplications, low-mappability regions and homopolymers. |
| Results affected | The explanation in report 01; no figure in it is wrong |
| Immediate action | Report 01 left as written; report 03 states the correction ("What changes from report 01") |
| Cause | Report 01 did not stratify accuracy by difficult region, and had one genome at one depth |
| Corrective action | Report 03 measures accuracy inside and outside the GIAB difficult-region strata (`08_strata.py`) |
| Effectiveness | Open: stratified accuracy is not yet a required step of a validation in SOP-01 section 8 |
| Risk register | R06, R14 |
| Status | open |

## NC-03: Pipeline runs lost to memory exhaustion

| | |
|---|---|
| Found | 2026-10-24 about 08:00 and 2026-10-25 06:49, report 03 runs |
| What | First, the machine ran out of memory while HG003 was aligning and three BAM to FASTQ conversions ran at the same time; every job stopped. Second, the HG005 alignment was terminated after 3 h 36 min during memory pressure. |
| Results affected | None. Interrupted tasks were rerun from their inputs and no partial file was used (`docs/reports/03-more-reference-samples/run_record.md`) |
| Immediate action | Partial FASTQ and temporary files deleted; runs restarted |
| Cause | Several heavy jobs at once on a 16 GB machine; the Nextflow memory limit covered only Nextflow's own tasks |
| Corrective action | One heavy job at a time; each process declares memory close to its real peak; executor limited to 9 GB; a memory monitor runs with long jobs; runs detached from the login session (design-decisions, "Running on a 16 GB laptop") |
| Effectiveness | HG004 and HG005 and every later run (reports 03 and 04, test data runs) completed without a memory failure |
| Risk register | R12 |
| Status | closed 2026-10-26 |

## NC-04: VCFEVAL can overwrite its truth input

| | |
|---|---|
| Found | 2026-10-26, while building the test data (PR #28) |
| What | The VCFEVAL process writes its working copy of the truth set to `truth.vcf.gz` in the task directory. A truth file supplied with that same name is staged as a link with that name, so the copy is written through the link and overwrites the input. The working copy of the calls, `calls.vcf.gz`, has the same weakness. |
| Results affected | None. The truth and call files used in reports 01 to 04 have other names. |
| Immediate action | Test truth files named `HG001_test_truth.*`; the fault is described in `tests/README.md` |
| Cause | Fixed working file names in the process script that can clash with staged input names |
| Corrective action | Planned: give the working copies names that cannot clash with an input. This edits workflow logic, so it goes through SOP-01 section 9; the test data check should show identical results |
| Effectiveness | Not yet checked |
| Risk register | R13 |
| Status | open |
