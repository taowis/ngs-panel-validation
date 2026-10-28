# Report 05 on NCI Gadi

DeepVariant has only x86-64 images, so report 05 ran it on NCI Gadi (project iz89) rather than on the arm64 Mac used for everything else. Everything that ran there is in this folder; nothing was written on Gadi by hand. The work folder was `/scratch/iz89/ktl590/ngs-panel-validation` (`$BASE`), which Gadi purges after 100 days without access.

1. Locally, cut the HG003 BAM to the SureSelect capture and copy it, the capture BED, the reference index and these scripts to `$BASE`: commands in `sync.sh` (`push`).
2. On Gadi, in `$BASE`: `qsub scripts/gadi/fetch.pbs` (copyq, which has internet). It pulls `docker://google/deepvariant:1.9.0` and `:1.10.0` to Singularity images and downloads hs37d5, checking it against the md5s in `data/MANIFEST.tsv`.
3. On Gadi: one `deepvariant.pbs` job per version and run (16 CPUs, normal queue), e.g. `qsub -v VERSION=1.9.0,RUN=run1,REGIONS=inputs/capture.HG003.bed scripts/gadi/deepvariant.pbs`. The jobs submitted, with their IDs and service units, are in `docs/reports/05-deepvariant-update/run_record.md`.
4. Locally, `sync.sh pull` copies the VCFs, logs and PBS output to `work_dv/gadi/`; `scripts/12_external_update.py` does the evaluation.
