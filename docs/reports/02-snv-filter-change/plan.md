# Change 02: SNV hard filter, plan

Written 2026-10-23, before any of the candidate filters below were run on any data.

## Reason for the change

The initial validation (`docs/reports/01-initial-validation/`) failed SNV sensitivity. The three panel SNVs missed in NIST7035 (5:80168937, 17:41245237 in BRCA1, 17:59763465 in BRIP1) were all called with the correct genotype and then removed by the SNP hard filter, each because SOR > 3. Exome-wide, about half the missed SNVs were called and filtered, mostly by the SOR rule alone and mostly homozygous ALT calls. With no reference reads at a homozygous site, SOR is driven by the strand balance of the ALT reads only, which in this capture is often one-sided near target edges.

These figures come from the whole exome, both chromosome halves below included. So the hypothesis (the SOR rule removes true calls) is not independent of the test set. What the split below keeps independent is the choice between candidates and the decision to accept the change.

## Candidates

The indel filter is not touched. Only the SNP_HARD expression in `GATHER_FILTER` changes.

- v1 (current): `QD<2 || FS>60 || MQ<40 || MQRankSum<-12.5 || ReadPosRankSum<-8 || SOR>3`
- A: as v1, but the SOR rule applies only to calls that are not homozygous ALT (`SOR>3 && GT!="AA"`)
- B: as v1, without the SOR rule

No other candidates will be tried on the development set. If neither is selected, the change is dropped and v1 stays.

## Data split

Regions are the exome evaluable region (Nextera capture & GIAB high-confidence) with the panel targets (`data/panels/panel_targets.bed`) removed.

- Development set: odd autosomes, NIST7035_run1 only. Used once, to choose between A and B.
- Test set: even autosomes, all three runs. Used once, after the choice is made and the pipeline is changed.
- Panel region: reported for v1 and v2, but its misses motivated the change, so it is not evidence that the change works. It is used for the review of discordant calls.

## Selection rule (development set)

Choose the candidate with the higher SNV F-measure (genotype mode, vcfeval), provided its SNV precision is no more than 0.005 below v1. If both fail the precision condition, select neither.

## Acceptance of the change (test set and full run)

The selected candidate is made the pipeline default (v2) and the full pipeline rerun on the same data into `results_v2/`. Following `config/acceptance.yaml`, change_control:

1. v1 itself failed several criteria, so "every criterion still passes" cannot apply as written. It is read as: no criterion that passed in v1 fails in v2.
2. On the test set, for each run, SNV and indel sensitivity and precision may not drop by more than 0.005 against v1. Indel metrics are expected to be identical.
3. Every call that differs between v1 and v2 inside the panel region is listed and reviewed against the truth set and the read evidence.
4. Within-run reproducibility (NIST7035 run twice) still gives identical VCF records.

The change is accepted only if 1 to 4 all hold. The outcome, accepted or rejected, is written up in this folder either way.
