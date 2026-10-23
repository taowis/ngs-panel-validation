"""Helpers shared by the report 03 scripts (06 to 10).

Each genome's inputs come from config/params/<genome>.yaml, the same file the
pipeline run used; its samples are the sample column of the samplesheet.
"""
import csv
import math
import os
import subprocess

import yaml

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
GENOMES = ["HG001", "HG002", "HG003", "HG004", "HG005"]
SDF = "data/reference/hs37d5.sdf"
REF = "data/reference/hs37d5.fa"
PANEL = "data/panels/panel_targets.bed"
PANEL_BY_GENE = "data/panels/panel_targets_by_gene.bed"
# HG001 has two libraries; report 03 uses NIST7035_run1 where one sample per genome is needed
MAIN_SAMPLE = {"HG001": "NIST7035_run1"}


def load(genome):
    g = yaml.safe_load(open(os.path.join(ROOT, f"config/params/{genome}.yaml")))
    with open(os.path.join(ROOT, g["samplesheet"])) as fh:
        g["samples"] = sorted({r["sample"] for r in csv.DictReader(fh)})
    g["main_sample"] = MAIN_SAMPLE.get(genome, genome)
    return g


def sh(cmd):
    return subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True, capture_output=True, text=True).stdout


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    ph = k / n
    d = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / d
    h = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)


def exome_bed(g, out):
    """capture & high-confidence, merged"""
    if not os.path.exists(out):
        sh(f"gzip -dc {g['capture']} | sed 's/^chr//' | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge "
           f"| bedtools intersect -a - -b {g['truth_bed']} | sort -k1,1 -k2,2n | bedtools merge > {out}")
    return out


def vcfeval(truth, calls, region, eval_regions, out, vtype, squash=False, threads=4):
    """SNVs or indels of both files, vcfeval over region; returns (tp_base, tp_call, fp, fn)."""
    sel = "snps" if vtype == "snv" else "indels"
    if not os.path.exists(f"{out}/summary.txt"):
        sh(f"rm -rf {out} {out}.*.vcf.gz* "
           f"&& bcftools view --no-version -v {sel} {truth} -Oz -o {out}.b.vcf.gz && tabix {out}.b.vcf.gz "
           f"&& bcftools view --no-version -v {sel} {calls} -Oz -o {out}.c.vcf.gz && tabix {out}.c.vcf.gz "
           f"&& rtg vcfeval -b {out}.b.vcf.gz -c {out}.c.vcf.gz -t {SDF} --bed-regions {region} "
           f"{'--evaluation-regions ' + eval_regions if eval_regions else ''} "
           f"{'--squash-ploidy' if squash else ''} --threads {threads} -o {out}")
    for line in open(f"{out}/summary.txt"):
        f = line.split()
        if f and f[0] == "None":
            return tuple(map(int, f[1:5]))
    return (0, 0, 0, 0)  # no variants in region


def rates(tpb, tpc, fp, fn):
    sens = tpb / (tpb + fn) if tpb + fn else float("nan")
    prec = tpc / (tpc + fp) if tpc + fp else float("nan")
    return sens, prec
