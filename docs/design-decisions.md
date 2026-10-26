# Design decisions

The main choices behind the pipeline and the validation, why each was made, and what the alternative was. The details and the evidence are in the reports and the SOP; this page points to them.

## Reference genome: GRCh37 (hs37d5)

The GIAB v4.2.1 benchmarks exist for both GRCh37 and GRCh38, but the capture files used here do not: the Nextera BED for HG001 and the SureSelect V5 baits published with the Oslo data are both in GRCh37 coordinates. Staying on GRCh37 avoids lifting capture regions over, which can drop or shift intervals at exactly the target edges where depth is already lowest. The cost is that GRCh37 lacks the GRCh38 fixes to some difficult regions. Report 03 measured accuracy in GIAB's difficult-region strata: segmental duplications are where most exome-wide SNVs are lost (sensitivity about 0.75 there, against 0.997 outside difficult regions).

## Truth sets and samples

HG001 (NA12878) came first because Garvan sequenced it as two separate libraries, which gives between-library reproducibility on real data, and the same FASTQ run twice gives within-run reproducibility. Its weakness is size: only 53 truth SNVs and no truth indels fall in the panel's evaluable region, so even a perfect result has a 95% lower bound of about 0.93.

Report 03 adds HG002 to HG005 for two reasons set out in its plan: more panel variants, and a test of report 01's explanation that the exome-wide failures came from the Garvan data (depth, capture) rather than the pipeline. The Oslo exomes use a different laboratory and capture kit, which is what that test needs. GIAB distributes them only as the laboratory's aligned BAMs, so the reads are taken back to FASTQ and run through the whole pipeline from alignment, rather than starting from Oslo's alignments. The laboratory's own calls on the same reads are then used as a comparison (report 03, analysis 8).

## Panel definition

Genes are the green (diagnostic grade) genes of the PanelApp Australia hereditary breast, ovarian and colorectal cancer panels, from snapshots kept in `data/panels/panelapp/`, so the gene list does not change under the analysis. Targets are the coding sequence of all protein-coding Ensembl 87 transcripts with 20 bp each side, merged. Taking every transcript rather than one canonical transcript per gene errs towards covering more: an exon that only some transcripts use is still checked for coverage. The 20 bp padding covers the splice sites next to each exon.

## Where accuracy is measured

Sensitivity and precision can only be measured where the truth set is confident and the capture kit has reads. So the primary region is panel targets & capture & GIAB high-confidence (`panel_evaluable.bed`, built per genome by `01_build_panel.py`). Panel bases outside it are not hidden: the capture gap (19% of the panel's coding sequence for the Nextera kit) is reported as a coverage failure.

The exome-wide region (capture & high-confidence) is a supporting region with its own criteria. It holds tens of thousands of variants, enough to see error patterns and to give tight intervals, which the 53 panel SNVs cannot.

## Acceptance criteria fixed in advance

`config/acceptance.yaml` was written on 2026-10-23 before any pipeline output existed, and any later change needs a written reason in a report. This is the point of a validation: the bar cannot be moved after the results are seen. Thresholds are set with reference to the AMP/CAP recommendations (Roy et al. 2018) but are this project's own. The point estimate must meet the threshold, and 95% Wilson intervals are always given next to it, because with few variants a pass on the point estimate can still leave a wide interval (SOP-01, section 8).

## Comparison engine: RTG vcfeval

vcfeval matches variants by haplotype, so the same variant written differently in the truth and the calls (left- or right-aligned indels, a complex change split into parts) still counts as a match. It is the engine that hap.py uses for GA4GH-style benchmarking; hap.py itself has no build for the arm64 machine this was run on. Each comparison runs in two modes: genotype mode (the genotype must match) and allele mode (`--squash-ploidy`, the allele only), so genotype errors can be told apart from missed variants. SNVs and indels are split with bcftools first so each has its own figures.

## Pipeline: BWA-MEM, single-sample HaplotypeCaller, hard filters

The pipeline is deliberately ordinary, so that the validation, not the tools, is the subject. Two choices matter for the results:

- BWA-MEM runs with `-K 100000000`, a fixed batch size, so the alignments do not depend on the number of threads. Together with pinned tool versions this is what makes the within-run check (same FASTQ run twice) give identical VCF records.
- Calls are filtered with GATK's generic hard filters, not VQSR. VQSR trains a model on the call set and needs more variants than a single exome provides; GATK's own guidance points small data sets to hard filters. The thresholds are GATK's published values, so the starting point was not tuned to these data.

Single-sample calling (no joint genotyping) matches how a diagnostic laboratory reports one patient at a time, and keeps each sample's result independent of what else was in the batch.

## Changing the pipeline: change 02 and the SOR rule

Report 01 failed and the diagnosis found one pipeline fault: the SOR rule of the SNP hard filter removed true homozygous ALT calls. The fix was handled as a controlled change, not an edit (`docs/reports/02-snv-filter-change/plan.md`):

- the plan, the two candidates and the selection rule were written before any candidate was run;
- candidates were compared on odd autosomes of one run (development set) and the chosen one tested once on even autosomes of all runs (test set), outside the panel, so the data used to decide was not the data used to judge;
- acceptance needed no criterion to fail that had passed before, no metric to drop by more than 0.005, every changed call in the panel to be reviewed, and within-run reproducibility to hold.

The panel misses that prompted the change were not used as evidence that it worked, since they were the reason for it. The HG002 to HG005 runs in report 03 played no part in the change and gave an independent check of it: v2 raised exome SNV sensitivity in all four genomes with at most 0.0003 lower precision, and recovered four panel SNVs that the SOR rule would have removed.

## Depth

The minimum of 20x over 98% of panel bases is a common laboratory figure, not one this project derived. Report 03 tests it by downsampling the HG002 BAM to 20x, 30x, 50x and 75x mean depth (fixed seed, so a rerun keeps the same reads) and running each level through the same calling and evaluation. The criterion is not changed inside that report. The result: panel SNVs were all found from 30x mean depth, exome SNV sensitivity stopped rising at about 75x, and even at 140x the SureSelect kit left 2.4% of panel bases below 20x, so meeting the 98% figure depends on the kit's evenness as well as on depth.

## Running on a 16 GB laptop

Nextflow with a pinned conda environment (`environment.yml`, full lock in `environment.lock.yml`) runs the same way on a laptop or a server; only the `local` profile in `nextflow.config` is machine-specific. The executor is limited to 9 GB, leaving room for the operating system, and each process declares memory close to its real peak (BWA-MEM about 7 GB with the index and sort buffers). After a run was lost on 2026-10-24 when several heavy jobs ran together, the report 03 genomes were run one heavy job at a time, with each genome's BAM to FASTQ conversion starting only after the previous genome had finished alignment. Interruptions and restarts are recorded as deviations in the run record.

## Containers and the test data

The validation runs used the conda environment on macOS. The `docker` and `apptainer` profiles run the same tool versions from a linux-64 image, for servers and clusters where conda is not wanted. The image installs from an explicit lock (`containers/linux-64.lock`, every package with its URL and md5), so rebuilding it later gives the same packages rather than whatever the solver picks that day. The profiles only change where the tools come from; `workflow/main.nf` and its commands are unchanged, so this is not a pipeline change under SOP-01 section 9.

A different build of the same versions could still behave differently, so the test data in `tests/` is the check. It is cut from the HG001 data (1.4 Mb of reference, 60,000 read pairs per library) so that GitHub Actions can run the whole pipeline in a few minutes, and its expected results come from the conda environment. CI requires the Docker and Apptainer runs to give the same VCF records and vcfeval counts. The test checks that nothing has changed; it says nothing about accuracy, which rests on the full-size runs in the reports.

## Out of scope

Copy number and exon-level deletions (GIAB has no exome-scale truth for them), structural variants, low allele fraction and mosaic variants, variant classification and clinical reporting (SOP-01, section 2). A laboratory using this pipeline would need separate validation for each of these.

## Review

SOP-01 asks for a second person to review validation reports and change decisions. In this project the author fills both roles, and the reports say so rather than imply an independent review.
