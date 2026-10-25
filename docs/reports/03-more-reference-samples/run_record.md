# Report 03: run record

Pipeline runs for HG002 to HG005, 2026-10-24 to 2026-10-25, on the 10-core, 16 GB Apple M4 used for reports 01 and 02, local executor. Each genome was run with the commands in `scripts/run_genomes.sh`, one genome at a time, into `results_v2_<genome>/`. The first line of each `pipeline_info/run.log` is the commit the run used:

| Genome | Commit | Pipeline code |
|---|---|---|
| HG003, HG002 | cb514bc | v2 as in c01b3a3 |
| HG005, HG004 | a14ed93 | same; the commits between cb514bc and a14ed93 change documentation and the project name only |

`workflow/main.nf` has since gained comments (PR #19) but no change to any command.

The Oslo BAMs were converted back to FASTQ (`scripts/00_download.sh oslo`) on the same machine. After the first memory failure, conversions and pipeline runs were scheduled so that only one heavy job ran at a time, with each genome's conversion starting after the previous genome had finished alignment, and a memory monitor running throughout.

## Timeline

| Time | Event |
|---|---|
| 2026-10-24 07:11 | HG003 started; HG002, HG004 and HG005 conversions running at the same time |
| about 08:00 | Deviation 1: the machine ran out of memory and every job stopped (HG003 at BWA_MEM, three conversions partial) |
| 15:25 | Partial FASTQ and temporary files deleted; restarted with one heavy job at a time |
| 15:26 | HG003 started again |
| 16:31 to 17:32 | HG002 conversion |
| 18:50 | HG003 done |
| 18:51 | HG002 started, including the four downsampled levels |
| 22:27 to 23:05 | HG005 conversion |
| 2026-10-25 01:38 | Run queue paused before HG005 and HG004 for a rename of the project folder; HG002 left running |
| 03:11 | HG002 done |
| 03:12 | Queue restarted for HG005 then HG004 |
| 06:49 | Deviation 2: the HG005 BWA_MEM task, Nextflow and the queue were terminated (SIGTERM, exit 143) after 3 h 36 min during a period of memory pressure; HG004 not started; HG005 FASTQ intact |
| 07:05 | Queue restarted for HG005 then HG004, detached from the login session |
| 08:08 to 08:44 | HG004 conversion |
| 09:41 | HG005 done |
| 12:35 | HG004 done |

All four runs finished with no failed tasks (HG002 178 tasks including the depth levels, the others 37 each). The partial work directory of the killed HG005 task was left in Nextflow's work folder and not reused.

## Effect on results

None expected. In both deviations the interrupted tasks were rerun from their inputs; no partial file was used. The rename moved the folder after HG002's outputs had been written, and the later runs used the new path throughout.

## Analyses

Run 2026-10-25 23:43 to 2026-10-26 01:49, one vcfeval at a time, with the memory monitor running and no warnings: `02_report.py`, `03_diagnose.py` for each genome, `06_pooled.py` to `10_interlab.py`, then `04_split_eval.py --set all` for each genome (analysis 3). `06_pooled.py` failed at the first attempt because `bcftools query` has no variant-type option; it was fixed (PR #24) and rerun. The `08_strata.py` read counts over PMS2 and MLH1 scan each whole BAM, which made that step take about 1.5 hours.

The exit codes the analysis queue wrote to its log were not reliable (they were taken after another command had run), so every output was checked by content instead: no error in the logs apart from the first `06_pooled.py` attempt; every table has the expected number of rows (for example 60 strata rows for five genomes, six strata and two variant types); the false-negative classes of each diagnosis sum to the allele-mode false negatives in that genome's vcfeval summary; the v2 rows of analysis 3 give the same counts as the pipeline's own vcfeval in `accuracy.tsv`; the pooled inputs equal the panel rows of `accuracy.tsv` and of report 02; and the strata, identity, sex, trio, inter-laboratory and analysis 3 tables match the rows the scripts printed as they ran.
