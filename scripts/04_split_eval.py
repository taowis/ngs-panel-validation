"""Accuracy on one half of the autosomes, for the filter change in
docs/reports/02-snv-filter-change/plan.md.

Region: exome evaluable (capture & GIAB high-confidence) minus the panel
targets, odd autosomes (dev) or even autosomes (test).

Two modes:
  --candidates v1 A B   refilter the v1 calls in <results>/vcf with each SNP
                        filter below and evaluate them
  --dirs results results_v2
                        evaluate the pipeline outputs in each directory as is

    python scripts/04_split_eval.py --set dev --candidates v1 A B --samples NIST7035_run1
    python scripts/04_split_eval.py --set test --dirs results results_v2
"""
import argparse
import csv
import os
import subprocess

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
p = argparse.ArgumentParser()
p.add_argument("--set", choices=["dev", "test"], required=True)
p.add_argument("--candidates", nargs="+")
p.add_argument("--dirs", nargs="+")
p.add_argument("--results", default="results")
p.add_argument("--samples", nargs="+", default=["NIST7035_run1", "NIST7035_run2", "NIST7086_run1"])
p.add_argument("--out", default="docs/reports/02-snv-filter-change")
args = p.parse_args()
os.chdir(ROOT)

SDF = "data/reference/hs37d5.sdf"
TRUTH = "data/truth/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED = "data/truth/HG001_GRCh37_1_22_v4.2.1_benchmark.bed"
CAPTURE = "data/raw/nexterarapidcapture_expandedexome_targetedregions.bed.gz"
PANEL = "data/panels/panel_targets.bed"

# must match GATHER_FILTER in workflow/main.nf (v1) and the plan (A, B)
SNP_RULES = {
    "v1": "INFO/QD<2 || INFO/FS>60 || INFO/MQ<40 || INFO/MQRankSum<-12.5 || INFO/ReadPosRankSum<-8 || INFO/SOR>3",
    "A": "INFO/QD<2 || INFO/FS>60 || INFO/MQ<40 || INFO/MQRankSum<-12.5 || INFO/ReadPosRankSum<-8 "
         "|| (INFO/SOR>3 && GT!=\"AA\")",
    "B": "INFO/QD<2 || INFO/FS>60 || INFO/MQ<40 || INFO/MQRankSum<-12.5 || INFO/ReadPosRankSum<-8",
}
INDEL_RULE = "INFO/QD<2 || INFO/FS>200 || INFO/ReadPosRankSum<-20 || INFO/SOR>10"


def sh(cmd):
    return subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True, capture_output=True, text=True).stdout


work = f"{args.results}/split_eval/{args.set}"
os.makedirs(work, exist_ok=True)
os.makedirs(args.out, exist_ok=True)
chroms = [str(c) for c in range(1 if args.set == "dev" else 2, 23, 2)]
region = f"{work}/region.bed"
keep = " || ".join(f'$1=="{c}"' for c in chroms)
sh(f"gzip -dc {CAPTURE} | sed 's/^chr//' | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge "
   f"| bedtools intersect -a - -b {TRUTH_BED} | bedtools subtract -a - -b {PANEL} "
   f"| awk '{keep}' | sort -k1,1 -k2,2n | bedtools merge > {region}")
for t, sel in [("snv", "snps"), ("indel", "indels")]:
    if not os.path.exists(f"{work}/truth.{t}.vcf.gz.tbi"):
        sh(f"bcftools view --no-version -v {sel} -R {region} {TRUTH} -Oz -o {work}/truth.{t}.vcf.gz "
           f"&& tabix -f {work}/truth.{t}.vcf.gz")

jobs = []  # (label, sample, vcf)
if args.candidates:
    for c in args.candidates:
        for s in args.samples:
            vcf = f"{work}/{c}.{s}.vcf.gz"
            sh(f"bcftools annotate --no-version -x FILTER {args.results}/vcf/{s}.vcf.gz "
               f"| bcftools filter --no-version -m+ -s SNP_HARD -e 'TYPE=\"snp\" && ({SNP_RULES[c]})' "
               f"| bcftools filter --no-version -m+ -s INDEL_HARD -e 'TYPE!=\"snp\" && ({INDEL_RULE})' "
               f"-Oz -o {vcf} && tabix -f {vcf}")
            if c == "v1":  # the refiltered v1 must reproduce the pipeline's own filtering
                a = sh(f"bcftools query -f '%CHROM %POS %REF %ALT %FILTER\\n' {args.results}/vcf/{s}.vcf.gz | md5sum")
                b = sh(f"bcftools query -f '%CHROM %POS %REF %ALT %FILTER\\n' {vcf} | md5sum")
                assert a == b, f"refiltered v1 differs from pipeline output for {s}"
            jobs.append((c, s, vcf))
else:
    for d in args.dirs:
        for s in args.samples:
            jobs.append((d, s, f"{d}/vcf/{s}.vcf.gz"))

rows = []
for label, s, vcf in jobs:
    for t, sel in [("snv", "snps"), ("indel", "indels")]:
        o = f"{work}/{os.path.basename(label)}.{s}.{t}"
        sh(f"rm -rf {o} && bcftools view --no-version -v {sel} {vcf} -Oz -o {o}.calls.vcf.gz "
           f"&& tabix -f {o}.calls.vcf.gz && rtg vcfeval -b {work}/truth.{t}.vcf.gz -c {o}.calls.vcf.gz -t {SDF} "
           f"--bed-regions {region} --evaluation-regions {TRUTH_BED} --threads 4 -o {o}")
        for line in open(f"{o}/summary.txt"):
            f = line.split()
            if f and f[0] == "None":
                tpb, tpc, fp, fn = map(int, f[1:5])
        sens, prec = tpb / (tpb + fn), tpc / (tpc + fp)
        rows.append([args.set, label, s, t, tpb, fn, tpc, fp, f"{sens:.4f}", f"{prec:.4f}",
                     f"{2 * sens * prec / (sens + prec):.4f}"])

name = f"{args.out}/{args.set}_{'candidates' if args.candidates else 'pipeline'}.tsv"
with open(name, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["set", "label", "sample", "type", "tp_baseline", "fn", "tp_call", "fp", "sensitivity", "precision",
                "f_measure"])
    w.writerows(rows)
for r in rows:
    print("\t".join(map(str, r)))
