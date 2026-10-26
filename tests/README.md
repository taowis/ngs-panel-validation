# Test data and CI

A small data set that runs the whole pipeline in a couple of minutes, used by the GitHub Actions workflow (`.github/workflows/test.yml`) and for checking an installation. It tests that the pipeline runs and gives the same calls as before. It is not a validation: the figures it produces are not judged against `config/acceptance.yaml`.

## What is in `tests/data/`

Built by `tests/make_test_data.py` from the full HG001 inputs and the pipeline v2 BAMs; rerunning it gives the same files.

- `ref/hs37d5.fa`: 24 contigs named 1 to 22, X and Y, as the pipeline expects, each one window of hs37d5 (1.4 Mb in all). Where a chromosome has a panel gene the window is that gene with 5 kb either side (MSH2 and EPCAM together on chromosome 2); otherwise it is 30 kb around the capture interval nearest the middle of the chromosome. `windows.tsv` gives each window's hs37d5 coordinates. The file keeps the name `hs37d5.fa` because the pipeline's commands use it.
- `NIST7035_R1/R2.fastq.gz`, `NIST7086_R1/R2.fastq.gz`: the read pairs from the two Garvan HG001 libraries with both mates in a window, about 60,000 pairs each, duplicates kept.
- `capture.bed.gz`, `panel_targets.bed`, `panel_evaluable.bed`, `HG001_test_truth.vcf.gz` and `.bed`: the Nextera capture, the panel BEDs and the NIST v4.2.1 HG001 benchmark, cut to the windows, with coordinates shifted to them.
- `samplesheet.csv`: NIST7035 twice (within-run check) and NIST7086 (between-library check).

The source data are GIAB's public HG001 files listed in `data/README.md`.

## Running it

The `test` profile points the pipeline at these files and lowers the resource requests; combine it with `local`, `docker` or `apptainer`. First index the reference (`tests/prepare_ref.sh`, inside the conda environment or the container), then:

```
nextflow run workflow/main.nf -profile local,test        # conda environment
nextflow run workflow/main.nf -profile docker,test       # image from nextflow.config
python3 tests/check_test_run.py                          # compare with tests/expected/
```

To use a locally built image: `docker build -t ngs-panel-validation:local .` and add `--container ngs-panel-validation:local`.

## What is checked

`tests/check_test_run.py` compares the run with `tests/expected/test_run.tsv`: no failed task, the VCF records of each sample (md5 of the lines after the header), the vcfeval counts for each sample, region, variant type and mode, and the within-run comparison. The expected values come from the conda environment on macOS (arm64). CI runs the same data in the linux-64 image under Docker and under Apptainer, so a pass also shows that the containers give the same calls as the environment used for reports 01 to 04.

A change that is meant to alter the calls will fail this test. Under SOP-01 section 9 that change has its own plan and report; once it is accepted, `python3 tests/check_test_run.py --write` updates the expected values, in the same pull request, with the reason given there.

## Found while building it

The VCFEVAL process writes its working copy of the truth set to `truth.vcf.gz` in the task directory. When the truth file supplied to it had that same name, the copy was written through the staged link and overwrote the input. The real truth files have other names, so no earlier run was affected; the test files are named `HG001_test_truth.*` for this reason.
