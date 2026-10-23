# Report NN: <title>

Pipeline version: commit <hash> (<main tool versions>). Run on <date> on <machine>, <executor>. Run logs: `<outdir>/pipeline_info/run.log`. Criteria: `config/acceptance.yaml`. Plan: `plan.md` (merged <date>, PR #n).

Data: <samples, libraries, sequencing, truth set>.

Deviations from the plan: <none, or what and why>.

## Outcome

One paragraph: passed or failed, and the criteria behind it.

| Criterion | <sample> | <sample> | Threshold | Result |
|---|---|---|---|---|

Intervals are 95% Wilson intervals.

## Why it fails

Only if it does. Separate pipeline causes from data causes, with the tables that show it.

## Limits

What this validation cannot show, and what data would be needed.

## Follow-up

Changes raised, further samples, open questions.

Reviewed by: <name, date>
