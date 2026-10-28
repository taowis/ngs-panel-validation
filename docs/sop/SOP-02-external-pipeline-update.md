# SOP-02: Verifying an update to an externally supplied analysis pipeline

| | |
|---|---|
| Document | SOP-02 |
| Version | 1.0 |
| Effective | 2026-10-28 |
| Author | Kaitao Lai |
| Next review | after its first use on a supplier's update, or 12 months |

This procedure was written for this project and first used in report 05. It is not the procedure of any accredited laboratory.

## 1. Purpose

To decide, before patient samples are run on it, whether a new version of an analysis pipeline that the laboratory did not write still meets the laboratory's acceptance criteria. The pipeline may be a commercial platform (for example SOPHiA DDM or Illumina DRAGEN), a supplier's cloud service, or a published caller such as DeepVariant. The laboratory does not see inside it, so the procedure works only from what the pipeline exports: VCF files, and BAM files where it produces them.

The supplier's own validation does not replace this. The laboratory has to show that the method performs as required on its own assay and data (ISO 15189:2022 clauses 6.8 and 7.3.2; NPAAC S6.8), and has to do so again for every version that can change results.

## 2. Scope

Any change made by a supplier that can alter calls or QC figures: a new software release, a new or retrained model, an updated reference, annotation or population database, changed default parameters or filters, or a change to how the supplier hosts the service. Changes the laboratory makes itself to its own pipeline follow SOP-01 section 9.

## 3. Assess the update

1. Get the release notes and record the old and new version identifiers exactly (version number, image digest or build number).
2. List what the notes say has changed and decide which outputs each change can affect: SNVs, indels, CNVs, filters, genotypes, annotation and HGVS names, QC metrics. Record the assessment in the change plan (`templates/change-plan.md`).
3. Check the update against `docs/risk-register.md`: which risks does it touch, and does it add one.
4. Decide the extent of verification from steps 2 and 3. A change to the caller or model needs at least the comparison in section 4. A change only to annotation needs the report-level comparison in section 5. If the notes are unclear, treat the update as a change to the caller.
5. Ask the supplier which data the new version was trained or tuned on. Reference samples used in training must not be used to judge it, or the result will be optimistic.

## 4. Compare the two versions on the same input

1. Use the same input for both versions: the same BAM or FASTQ, the same capture regions, the same reference. Hold everything else fixed.
2. Use reference samples with a truth set (GIAB) that the new version was not trained on.
3. Run each version twice on the same input. The repeat shows how much results move between runs of the same version; a difference between versions smaller than that is not evidence of a change.
4. Export the VCFs and run `scripts/12_external_update.py` with the old version as `--before` and the new one as `--after`. It scores each version against the truth in the panel and exome evaluable regions, applies the change-control rules in `config/acceptance.yaml`, and lists every PASS call in the panel targets that differs between the versions.
5. Review every differing panel call against the truth and the read evidence, as in SOP-01 section 9.6.

## 5. Report-level comparison

Where the platform also annotates or classifies variants, compare the reported variants on the same samples: gene, transcript, HGVS c. and p. descriptions, zygosity and classification. Any difference is reviewed and explained before release.

## 6. Decision

The update is accepted only if every rule in section 4 passes, every difference has been reviewed, and section 5 shows no unexplained change. Record the decision in a report from `templates/change-report.md`, with the supplier's version identifiers, the job or run records, and the cost where the runs used shared computing.

If the update fails, it is not used for patient samples. Raise it with the supplier with a report they can reproduce (input, versions, the calls that changed and why they are wrong), and enter it in `docs/nonconformities.md`.

## 7. Records

Release notes, the change plan, the VCFs of every run (or their md5 and location), the outputs of `12_external_update.py` and the report, kept as in SOP-01 section 11.

## 8. References

International Organization for Standardization. ISO 15189:2022 Medical laboratories - Requirements for quality and competence. 4th ed. Geneva: ISO; 2022.

National Pathology Accreditation Advisory Council. Requirements for Human Medical Genome Testing Utilising Massively Parallel Sequencing Technologies (First Edition 2017). Canberra: Commonwealth of Australia; 2017.
