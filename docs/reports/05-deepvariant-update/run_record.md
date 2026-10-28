# Report 05: run record

NCI Gadi, project iz89, work folder `/scratch/iz89/ktl590/ngs-panel-validation`. Scripts in `scripts/gadi/`, copied to Gadi with `scripts/gadi/sync.sh push`; nothing was edited on Gadi.

## Inputs

| Input | Detail |
|---|---|
| BAM | `results_v2_HG003/bam/HG003.md.bam` (pipeline v2, report 03) cut to the merged SureSelect V5 capture with `samtools view -M -L`: 6,456,911,851 bytes, md5 `eb143f82979327c136c077015d537d96`, the same on the Mac and on Gadi after transfer |
| Regions | `capture.HG003.bed`: 194,819 intervals, 70.6 Mb |
| Reference | hs37d5 downloaded on Gadi; `hs37d5.fa.gz` and `hs37d5.fa` md5 checked against `data/MANIFEST.tsv` (`gadi/dv_fetch.o182417116`) |
| Images | Singularity images pulled from Docker Hub; md5 in `gadi/images.md5`, digests in `data/MANIFEST.tsv` |

## Jobs

All exited 0. DeepVariant jobs: normal queue, 16 CPUs, 64 GB requested, `run_deepvariant --model_type=WES --num_shards=16`.

| Job | What | Node | Walltime | Memory used | SU |
|---|---|---|---|---|---|
| 182417116 | fetch (copyq) | | 00:20:56 | | 2.79 |
| 182419601 | 1.9.0 run1 | gadi-cpu-clx-2218 | 00:18:19 | 19.2 GB | 9.77 |
| 182421035 | 1.9.0 run2 | gadi-cpu-clx-1890 | 00:18:04 | 19.3 GB | 9.64 |
| 182421034 | 1.10.0 run1 | gadi-cpu-clx-2587 | 00:22:01 | 22.4 GB | 11.74 |
| 182421036 | 1.10.0 run2 | gadi-cpu-clx-1122 | 00:19:49 | 22.3 GB | 10.57 |

Total 44.5 SU. The iz89 balance read 3.27 KSU used before the fetch job and 3.38 KSU after the last run (the balance includes other users' jobs).

## Outputs

VCFs are not committed; they are in `work_dv/gadi/results/` locally and on Gadi until the scratch purge, and can be made again with the scripts above.

| Run | `HG003.vcf.gz` md5 | Bytes | md5 of the records (header excluded) |
|---|---|---|---|
| 1.9.0 run1 | 7bdf6c5b05d91f369e4d3ae6210a2e6d | 1,650,506 | 1149c1d911c087040e37bcc1bca08108 |
| 1.9.0 run2 | 7bdf6c5b05d91f369e4d3ae6210a2e6d | 1,650,506 | 1149c1d911c087040e37bcc1bca08108 |
| 1.10.0 run1 | 2ebb7dfe3bddcbb59499fb1a1edf9f8f | 1,651,394 | f2a46bcd045a89b96618bfd83dd85e76 |
| 1.10.0 run2 | 2ebb7dfe3bddcbb59499fb1a1edf9f8f | 1,651,394 | f2a46bcd045a89b96618bfd83dd85e76 |

## Evaluation

On the Mac, in the conda environment:

```
R=work_dv/gadi/results
python scripts/12_external_update.py --genome HG003 \
    --before 1.9.0=$R/deepvariant_1.9.0_run1/HG003.vcf.gz,$R/deepvariant_1.9.0_run2/HG003.vcf.gz \
    --after 1.10.0=$R/deepvariant_1.10.0_run1/HG003.vcf.gz,$R/deepvariant_1.10.0_run2/HG003.vcf.gz \
    --context gatk_v2=results_v2_HG003/vcf/HG003.vcf.gz \
    --out docs/reports/05-deepvariant-update
```

The WES model files were compared inside the two images on a Gadi login node (`singularity exec … md5sum` over `/opt/models/wes`, output in `gadi/wes_model_md5.txt`).

## Deviation

DeepVariant was first pulled with Docker on the Mac. The 1.9.0 image unpacked to 10.3 GB, leaving 1.8 GB free on the system disk while swap stood at 14 GB, and the image is x86-64 only, so it would have run under emulation. It was removed without being run and the work moved to Gadi. No result came from the Mac attempt.
