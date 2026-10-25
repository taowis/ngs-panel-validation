# Report 03: more reference samples

Pipeline: v2 as accepted in report 02, unchanged. Runs from commits cb514bc (HG003, HG002) and a14ed93 (HG005, HG004), which have the same pipeline code; dates, interruptions and restarts are in `run_record.md`. Criteria: `config/acceptance.yaml`, unchanged. The plan (`plan.md`) and its additions were written on 2026-10-24 before any of these samples were run.

Data: GIAB HG002, HG003 and HG004 (Ashkenazi trio) and HG005 (Han Chinese son), exomes sequenced at Oslo University Hospital with Agilent SureSelect V5 on a HiSeq 2500, 2 x 126 bp. GIAB distributes only the laboratory's aligned BAMs; the reads were taken back to FASTQ and run through the whole pipeline from alignment. Truth: NIST v4.2.1 GRCh37 (the `_noinconsistent` high-confidence BED for the trio). Mean depth over the capture is 140x (HG002), 120x (HG003), 132x (HG004) and 114x (HG005), against 47x for HG001 NIST7035.

Tables in this folder come from `02_report.py`, `03_diagnose.py` (one set per genome), `04_split_eval.py --set all` and `06_pooled.py` to `10_interlab.py`.

## Outcome

v2 does not meet every acceptance criterion on these genomes, and the reasons differ from those report 01 gave.

Panel SNV sensitivity passes in all four genomes, and over the five genomes pooled no panel SNV is missed (316/316). Panel SNV precision fails in HG002 and HG003 because of two adjacent false calls in PMS2. Exome-wide SNV and indel criteria fail in every genome, as they did for HG001, although the depth is now two to three times higher. Outside GIAB's difficult regions they pass: SNV sensitivity is 0.996 to 0.998. The exome-wide failures are therefore concentrated in segmental duplications, low-mappability regions and homopolymers, which report 01 did not consider. The sequencing laboratory's own calls on the same reads also fail the exome criteria, with lower SNV sensitivity (0.955) and higher SNV precision than this pipeline. Panel coverage misses 98% at 20x narrowly in every genome, mostly in PMS2, CHEK2 and STK11.

Change 02 holds up on data that played no part in it: removing the SOR rule raised exome SNV sensitivity by 4 to 5 points in every genome, cost at most 0.0003 in precision, and recovered four panel SNVs that v1 would have missed.

## Analysis 1: each genome against the criteria

`acceptance_checks.tsv`, `accuracy.tsv`. Intervals are 95% Wilson intervals.

| Criterion | HG002 | HG003 | HG004 | HG005 | Threshold |
|---|---|---|---|---|---|
| Panel SNV sensitivity | 63/63 (0.943-1) | 72/72 (0.949-1) | 59/59 (0.939-1) | 69/69 (0.947-1) | 0.99 |
| Panel SNV precision | 63/65 = 0.969 | 72/74 = 0.973 | 1.000 | 1.000 | 0.99 |
| Panel indels, truth variants | 4 | 3 | 3 | 1 | too few to judge |
| Exome SNV sensitivity | 0.981 | 0.982 | 0.981 | 0.979 | 0.99 |
| Exome SNV precision | 0.994 | 0.993 | 0.994 | 0.993 | 0.99 |
| Exome SNV genotype concordance | 0.999 | 0.999 | 0.999 | 0.998 | 0.99 |
| Exome indel sensitivity / precision | 0.899 / 0.904 | 0.894 / 0.895 | 0.900 / 0.898 | 0.938 / 0.940 | 0.95 |
| Exome indel genotype concordance | 0.961 | 0.955 | 0.958 | 0.974 | 0.99 |
| Panel bases at >= 20x | 0.976 | 0.976 | 0.975 | 0.974 | 0.98 |

Panel SNV sensitivity passes on the point estimate in all four, but with 59 to 72 truth SNVs per genome no single genome's lower bound reaches 0.99. The reproducibility criteria are not assessable: each genome has one library and one run.

### The panel false calls

Every false call in the panel evaluable region, from vcfeval in genotype mode:

| Position | Gene | Call | Genomes | Note |
|---|---|---|---|---|
| 7:6026707 T>C and 7:6026708 C>T | PMS2 | 0/1, PASS | HG002, HG003 | adjacent pair, about 19% ALT reads (79 of 416 in HG002), QD 5.4, MQ 57; not in the truth |
| 11:108114661 AT>A | ATM | 0/1 | HG003, HG004, HG005 | next to a run of 14 T; depth 14 to 18 |
| 11:108121410 C>CT | ATM | 0/1 | HG002, HG003 | truth has this insertion; next to a run of 15 T; depth 7 to 8 |
| 2:215657182 TAA>T | BARD1 | 1/1 (HG002), 1/2 (HG004) | HG002, HG004 | truth has TA>T; run of 16 A |

The PMS2 pair is the only reason panel SNV precision fails. A low allele fraction, two changes side by side and slightly reduced mapping quality are what reads from the pseudogene PMS2CL look like when placed on PMS2. The pair is seen in the father and the son but not the mother, so it may be a sequence difference carried by that family rather than random error; it cannot be settled without a PMS2-specific method such as long-range PCR. The indel false calls are all at homopolymers of 14 bp or more, at low depth, and at two of the three positions the truth has an indel there too, so the call has the wrong length or genotype rather than being invented.

## Why the exome-wide criteria fail

Report 01 put HG001's exome-wide failures down to shallow data, capture gaps and homopolymer indels. The Oslo data are much deeper and still fail, so depth cannot be the main reason. Four analyses together locate the cause.

Difficult regions (analysis 6, `strata_accuracy.tsv`): exome evaluable region split by GIAB v3.6 stratification.

| Stratum | Share of evaluable bp | SNV sensitivity | SNV precision | Indel sensitivity | Indel precision |
|---|---|---|---|---|---|
| Not in any difficult region | about 78% | 0.996-0.998 | 0.998-0.999 | 0.991-0.994 | 0.987-0.992 |
| All difficult regions | about 22% | 0.928-0.939 | 0.978-0.980 | 0.839-0.899 | 0.842-0.904 |
| Segmental duplications | 3.5% | 0.743-0.760 | 0.922-0.924 | 0.775-0.805 | 0.831-0.880 |
| Homopolymers | 1.5% | 0.969-0.975 | 0.964-0.978 | 0.796-0.912 | 0.790-0.893 |
| GC below 25% or above 65% | 16% | 0.971-0.980 | 0.990-0.991 | 0.915-0.939 | 0.930-0.956 |

Ranges are over HG002 to HG005. Outside the difficult regions the exome-wide sensitivity and precision criteria pass in every genome. For HG001 the same split gives SNV sensitivity of 0.926 outside the difficult regions, so for the Garvan data the shallow depth found in report 01 is real and adds to the difficult-region losses.

Diagnosis (`diagnosis_*.HG00x.tsv`): of HG002's 1,005 missed SNVs (allele level), 426 were called and filtered, mostly by MQ < 40 (the mapping-quality rule, which marks reads that could belong elsewhere in the genome), and 579 were not called. 459 of the misses are below 20x despite a mean of 140x, so the low depth is local, not general. 221 of HG002's 372 false indels are 1 bp deletions next to homopolymers of 8 bp or longer, the pattern report 01 found. The other genomes look the same.

Depth (analysis 5): see below. Exome SNV sensitivity stops rising at about 75x.

The sequencing laboratory (analysis 8): see below. Its own calls do not meet the exome criteria either.

## Analysis 2: panel accuracy pooled over five genomes

`pooled_panel_accuracy.tsv`, `pooled_panel_inputs.tsv`. HG001 is NIST7035_run1 from report 02.

| Type | TP | FN | FP | Sensitivity (95% CI) | Precision (95% CI) | Distinct truth variants |
|---|---|---|---|---|---|---|
| SNV | 316 | 0 | 4 | 1.000 (0.988-1.000) | 0.988 (0.968-0.995) | 116 |
| Indel | 7 | 4 | 7 | 0.64 (0.35-0.85) | 0.50 (0.27-0.73) | 5 |

The 316 truth SNVs are only 116 distinct variants, because the trio share alleles and common variants recur across genomes. The Wilson interval treats all 316 as independent and is therefore optimistic. Pooled SNV precision, 0.988, is just below 0.99, and all four false positives are the PMS2 pair. Five distinct indels are still too few to say anything about panel indel calling.

## Analysis 3: independent check of change 02

The v1 SNP filter (with SOR > 3) and v2 applied to the same calls (`all_candidates.HG00x.tsv`; `04_split_eval.py` checks that its v2 refiltering reproduces the pipeline's output exactly before it is used).

| Genome | Exome SNV sensitivity v1 → v2 | Exome SNV precision v1 → v2 | Panel SNV found v1 → v2 |
|---|---|---|---|
| HG002 | 0.9385 → 0.9814 | 0.9944 → 0.9941 | 63/63 → 63/63 |
| HG003 | 0.9372 → 0.9820 | 0.9934 → 0.9933 | 71/72 → 72/72 |
| HG004 | 0.9390 → 0.9809 | 0.9940 → 0.9938 | 58/59 → 59/59 |
| HG005 | 0.9293 → 0.9793 | 0.9933 → 0.9932 | 67/69 → 69/69 |

Indel results are identical, as expected. The plan's condition (v2 sensitivity at least v1's and precision no more than 0.005 lower in every genome) holds in all four; the largest precision drop is 0.0003, much smaller than on HG001. The four panel SNVs v1 would have missed are 17:33433487 (RAD51D), 7:6038722 (PMS2), 2:47693959 (MSH2) and 17:59763114 (BRIP1). Each is homozygous ALT with no reference reads, QD 28 to 34, MQ 60, and failed only SOR > 3, the same pattern as the three HG001 misses that prompted the change.

## Analysis 4: known pathogenic sites

`clinvar_sites_summary.tsv`, `clinvar_sites_by_gene.tsv`, `clinvar_sites_below_20x.tsv`. ClinVar 2026-10-19, GRCh37, pathogenic or likely pathogenic, not conflicting, up to 50 bp: 38,833 sites in the panel targets (820 longer records left out).

| Sample | Sites at >= 20x |
|---|---|
| HG001 NIST7035 / NIST7086 | 0.871 / 0.886 |
| HG002 to HG005 | 0.981-0.983 |

On the Oslo data the genes with the lowest fraction are PMS2 (0.77-0.81 of 884 sites), CHEK2 (0.88-0.89 of 1,005) and NF1 (0.93-0.94 of 5,488); all others are at 0.97 or above except EPCAM in HG005 (43 of 46). Every site below 20x is listed. Several hundred sites outside the capture BED still reach 20x from reads that run past the baits, which is why sites at 20x can outnumber sites in the capture.

## Analysis 5: depth

HG002 downsampled to a mean of 20x, 30x, 50x and 75x over the capture (`depth_downsampling.tsv`).

| Mean depth | Panel bases >= 20x | Exome SNV sensitivity / precision | Exome indel sensitivity / precision | Panel SNV found |
|---|---|---|---|---|
| 20x | 0.677 | 0.920 / 0.980 | 0.731 / 0.869 | 62/63 |
| 30x | 0.821 | 0.952 / 0.987 | 0.792 / 0.888 | 63/63 |
| 50x | 0.921 | 0.972 / 0.991 | 0.844 / 0.898 | 63/63 |
| 75x | 0.956 | 0.978 / 0.993 | 0.867 / 0.900 | 63/63 |
| 140x (full) | 0.976 | 0.981 / 0.994 | 0.899 / 0.904 | 63/63 |

Panel SNVs are all found from 30x. Exome SNV sensitivity gains only 0.003 between 75x and 140x, so the remaining misses are not a matter of depth. With this kit, even 140x leaves 2.4% of panel bases below 20x, which says the 98% criterion depends on how evenly a kit covers the panel as much as on how deeply it is sequenced. A mean depth of 20x is far from enough for the 20x-per-base criterion. As the plan says, the criterion is not changed here.

## Analysis 6: PMS2

`pms2.tsv`, with MLH1 for comparison.

| Gene | Coding bp in high-confidence | Fraction >= 20x | Reads with MAPQ < 20 |
|---|---|---|---|
| PMS2 | 0.95-0.97 | 0.77-0.79 | 13-14% |
| MLH1 | 0.985-0.986 | 0.995-0.996 | 0% |

Figures are for HG002 to HG005; HG001 has 17% of PMS2 reads at MAPQ < 20. About a fifth of PMS2's coding sequence stays below 20x and many of its reads map ambiguously because of PMS2CL, and the PMS2 false calls above come from the same source. This is a limitation of the method, not something a filter can fix; diagnostic laboratories test PMS2 with a pseudogene-aware method.

## Analysis 7: sample identity

`identity.tsv`, `sex_check.tsv`, `trio_mendelian.tsv`.

Each sample's PASS SNVs on chromosomes 20 to 22 agree best with its own genome's truth (0.987 to 0.993 for the Oslo samples, 0.94 for HG001) and at 0.33 to 0.42 with unrelated genomes. HG002 agrees at 0.55 to 0.60 with its parents, as a child should. Inferred sex matches for all six samples (chrX heterozygous fraction 0.59 to 0.63 and chrY depth 0.001 to 0.019 of the capture mean in females; 0.03 to 0.04 and 0.34 to 0.35 in males). In the trio, at positions at >= 20x in all three, 0.66% of SNV sites and 7.9% of indel sites are Mendelian inconsistent; an error in any one of the three counts, and the indel figure reflects the homopolymer errors above.

## Analysis 8: the sequencing laboratory's calls

`interlab_accuracy.tsv`, `interlab_concordance.tsv`. Oslo's GATK HaplotypeCaller calls with VQSR, published by GIAB (joint calls for the trio, single-sample for HG005), scored over the same exome evaluable region.

| Genome | SNV sensitivity, this / Oslo | SNV precision, this / Oslo | Indel sensitivity, this / Oslo | Indel precision, this / Oslo |
|---|---|---|---|---|
| HG002 | 0.981 / 0.955 | 0.994 / 0.999 | 0.899 / 0.893 | 0.904 / 0.900 |
| HG003 | 0.982 / 0.955 | 0.993 / 0.998 | 0.894 / 0.886 | 0.895 / 0.897 |
| HG004 | 0.981 / 0.954 | 0.994 / 0.998 | 0.900 / 0.890 | 0.898 / 0.899 |
| HG005 | 0.979 / 0.956 | 0.993 / 0.999 | 0.938 / 0.914 | 0.940 / 0.960 |

This pipeline finds 1,360 to 1,550 more true SNVs per genome and makes 250 to 310 more false SNV calls; VQSR filters harder. Indel results are close. The two call sets agree on 0.966 to 0.969 of SNV sites, and most disagreements are calls made here and not by Oslo (about 1,800 against 50). There is no gap pointing to a fault in this pipeline, and the laboratory's own calls also fail the exome criteria on these data. The Oslo calls were made in 2015 with the tools of that time, so this is a comparison with one laboratory's practice, not with a current standard.

## What changes from report 01

Report 01 explained HG001's exome-wide failures by depth, capture and homopolymer indels. On HG001 the depth part is still supported: even outside difficult regions HG001 reaches only 0.926 SNV sensitivity. But the Oslo data show that at diagnostic-like depth the exome-wide criteria still fail, and that the failures sit in segmental duplications, low-mappability regions and homopolymers. Report 01 did not look at these regions, and its statement that the remaining failures "need data at diagnostic depth" was incomplete. Report 01 is left as written; this report corrects it.

## Limits

- Panel indels remain untested in practice: five distinct truth indels over five genomes.
- Pooled panel figures count shared variants more than once; the intervals are optimistic.
- The Oslo exomes were sequenced in 2015, one lane each, with one capture kit. They are a second laboratory, not a current diagnostic panel.
- Exome-wide figures include regions a diagnostic panel would not report without other methods; the criteria in `acceptance.yaml` were written for the whole evaluable region and are applied as written.

## Decision and follow-up

Report 03 makes no change to the pipeline or the criteria. It finds:

1. Change 02 is confirmed on independent data.
2. Panel SNV sensitivity meets the criterion in every genome; panel SNV precision fails in two genomes because of a PMS2 pseudogene artefact.
3. The exome-wide criteria fail because of difficult regions, on this pipeline and on the laboratory's own calls, and pass outside them.
4. Panel coverage misses 98% at 20x narrowly, with PMS2, CHEK2 and STK11 the main gaps.

Each of the following would be a new change with its own plan:

- Judge exome-wide accuracy outside GIAB difficult regions and report the difficult strata separately, with a written reason in `acceptance.yaml`. This would make the criteria match what short-read exome data can show; it must not be presented as the pipeline improving.
- Exclude the PMS2 region covered by PMS2CL from short-read reporting, or flag it for a pseudogene-aware method, in SOP-01.
- Flag indels next to homopolymers of 8 bp or more for confirmation.
