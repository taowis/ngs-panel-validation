"""Accuracy on one half of the autosomes, for the filter change in
docs/reports/02-snv-filter-change/plan.md.

Region: exome evaluable (capture & GIAB high-confidence) minus the panel
targets, odd autosomes (dev) or even autosomes (test).

--set all is for genomes that took no part in choosing the change (report 03):
all autosomes of the exome evaluable region, panel included, and also the
panel evaluable region given by --panel-evaluable.

Two modes:
  --candidates v1 A B   refilter the v1 calls in <results>/vcf with each SNP
                        filter below and evaluate them
  --dirs results results_v2
                        evaluate the pipeline outputs in each directory as is

    python scripts/04_split_eval.py --set dev --candidates v1 A B --samples NIST7035_run1
    python scripts/04_split_eval.py --set test --dirs results results_v2
    python scripts/04_split_eval.py --set all --candidates v1 B --pipeline-filter B \
        --results results_v2_HG002 --samples HG002 --tag HG002 --truth ... --truth-bed ... \
        --capture ... --panel-evaluable data/panels/panel_evaluable.HG002.bed \
        --out docs/reports/03-more-reference-samples
"""
import argparse
import csv
import os
import subprocess

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
p = argparse.ArgumentParser()
p.add_argument("--set", choices=["dev", "test", "all"], required=True)
p.add_argument("--candidates", nargs="+")
p.add_argument("--dirs", nargs="+")
p.add_argument("--results", default="results")
p.add_argument("--samples", nargs="+", default=["NIST7035_run1", "NIST7035_run2", "NIST7086_run1"])
p.add_argument("--out", default="docs/reports/02-snv-filter-change")
p.add_argument("--truth", default="data/truth/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz")
p.add_argument("--truth-bed", default="data/truth/HG001_GRCh37_1_22_v4.2.1_benchmark.bed")
p.add_argument("--capture", default="data/raw/nexterarapidcapture_expandedexome_targetedregions.bed.gz")
p.add_argument("--panel-evaluable", help="with --set all")
p.add_argument("--pipeline-filter", default="v1", help="candidate that must reproduce the pipeline's filtering")
p.add_argument("--tag", help="suffix for the work folder and output table")
args = p.parse_args()
os.chdir(ROOT)

SDF = "data/reference/hs37d5.sdf"
TRUTH, TRUTH_BED, CAPTURE = args.truth, args.truth_bed, args.capture
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


sfx = f".{args.tag}" if args.tag else ""
work = f"{args.results}/split_eval/{args.set}{sfx}"
os.makedirs(work, exist_ok=True)
os.makedirs(args.out, exist_ok=True)
if args.set == "all":
    chroms = [str(c) for c in range(1, 23)]
    exclude = ""
else:
    chroms = [str(c) for c in range(1 if args.set == "dev" else 2, 23, 2)]
    exclude = f"| bedtools subtract -a - -b {PANEL} "
keep = " || ".join(f'$1=="{c}"' for c in chroms)
regions = {"exome" if args.set == "all" else args.set: f"{work}/region.bed"}
sh(f"gzip -dc {CAPTURE} | sed 's/^chr//' | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge "
   f"| bedtools intersect -a - -b {TRUTH_BED} {exclude}"
   f"| awk '{keep}' | sort -k1,1 -k2,2n | bedtools merge > {work}/region.bed")
if args.set == "all":
    regions["panel"] = args.panel_evaluable
for t, sel in [("snv", "snps"), ("indel", "indels")]:
    if not os.path.exists(f"{work}/truth.{t}.vcf.gz.tbi"):
        sh(f"bcftools view --no-version -v {sel} -R {work}/region.bed {TRUTH} -Oz -o {work}/truth.{t}.vcf.gz "
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
            if c == args.pipeline_filter:  # must reproduce the pipeline's own filtering
                a = sh(f"bcftools query -f '%CHROM %POS %REF %ALT %FILTER\\n' {args.results}/vcf/{s}.vcf.gz | md5sum")
                b = sh(f"bcftools query -f '%CHROM %POS %REF %ALT %FILTER\\n' {vcf} | md5sum")
                assert a == b, f"refiltered {c} differs from pipeline output for {s}"
            jobs.append((c, s, vcf))
else:
    for d in args.dirs:
        for s in args.samples:
            jobs.append((d, s, f"{d}/vcf/{s}.vcf.gz"))

rows = []
for label, s, vcf in jobs:
    for t, sel in [("snv", "snps"), ("indel", "indels")]:
        o = f"{work}/{os.path.basename(label)}.{s}.{t}"
        sh(f"bcftools view --no-version -v {sel} {vcf} -Oz -o {o}.calls.vcf.gz && tabix -f {o}.calls.vcf.gz")
        for rname, bed in regions.items():
            ro = f"{o}.{rname}"
            sh(f"rm -rf {ro} && rtg vcfeval -b {work}/truth.{t}.vcf.gz -c {o}.calls.vcf.gz -t {SDF} "
               f"--bed-regions {bed} --evaluation-regions {TRUTH_BED} --threads 4 -o {ro}")
            tpb = tpc = fp = fn = 0  # no variants in region
            for line in open(f"{ro}/summary.txt"):
                f = line.split()
                if f and f[0] == "None":
                    tpb, tpc, fp, fn = map(int, f[1:5])
            sens = tpb / (tpb + fn) if tpb + fn else float("nan")
            prec = tpc / (tpc + fp) if tpc + fp else float("nan")
            fm = 2 * sens * prec / (sens + prec) if sens + prec else float("nan")
            rows.append([rname, label, s, t, tpb, fn, tpc, fp, f"{sens:.4f}", f"{prec:.4f}", f"{fm:.4f}"])

name = f"{args.out}/{args.set}_{'candidates' if args.candidates else 'pipeline'}{sfx}.tsv"
with open(name, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["set", "label", "sample", "type", "tp_baseline", "fn", "tp_call", "fp", "sensitivity", "precision",
                "f_measure"])
    w.writerows(rows)
for r in rows:
    print("\t".join(map(str, r)))
