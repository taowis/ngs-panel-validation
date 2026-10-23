"""Explain the exome-wide misses and false calls of a pipeline run.

For each sample, reruns vcfeval in allele mode over the exome evaluable region
(capture & GIAB high-confidence) and classifies:
  false negatives  by whether a filtered call sits at the same position, which
                   hard-filter rule(s) removed it, its genotype, and whether
                   the position is callable (>= min_depth) in that sample
  false positives  (indels) by the longest homopolymer next to the call

Run inside the conda env (needs rtg, bcftools, bedtools, samtools):

    python scripts/03_diagnose.py --name 01-initial-validation

For another genome, also give --truth, --truth-bed and --capture.
"""
import argparse
import bisect
import csv
import os
import subprocess
from collections import Counter

import yaml

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
p = argparse.ArgumentParser()
p.add_argument("--name", required=True)
p.add_argument("--results", default="results")
p.add_argument("--samples", nargs="+", default=["NIST7035_run1", "NIST7086_run1"])
p.add_argument("--truth", default="data/truth/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz")
p.add_argument("--truth-bed", default="data/truth/HG001_GRCh37_1_22_v4.2.1_benchmark.bed")
p.add_argument("--capture", default="data/raw/nexterarapidcapture_expandedexome_targetedregions.bed.gz")
p.add_argument("--tag", help="suffix for the output tables, e.g. the genome")
args = p.parse_args()
os.chdir(ROOT)
cfg = yaml.safe_load(open("config/acceptance.yaml"))
out = f"docs/reports/{args.name}"
work = f"{args.results}/diagnose"
sfx = f".{args.tag}" if args.tag else ""
os.makedirs(out, exist_ok=True)
os.makedirs(work, exist_ok=True)

REF = "data/reference/hs37d5.fa"
SDF = "data/reference/hs37d5.sdf"
TRUTH, TRUTH_BED, CAPTURE = args.truth, args.truth_bed, args.capture

# the v1 SNP rules from GATHER_FILTER (workflow/main.nf), before change 02
SNP_RULES = {"QD": lambda v: v["QD"] < 2, "FS": lambda v: v["FS"] > 60, "MQ": lambda v: v["MQ"] < 40,
             "MQRankSum": lambda v: v["MQRankSum"] < -12.5, "ReadPosRankSum": lambda v: v["ReadPosRankSum"] < -8,
             "SOR": lambda v: v["SOR"] > 3}


def sh(cmd):
    return subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True, capture_output=True, text=True).stdout


def num(x):
    return float(x) if x not in (".", "") else float("nan")


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


exome = f"{work}/exome_evaluable.bed"
sh(f"gzip -dc {CAPTURE} | sed 's/^chr//' | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge "
   f"| bedtools intersect -a - -b {TRUTH_BED} | sort -k1,1 -k2,2n | bedtools merge > {exome}")

fn_rows, rule_rows, fp_rows = [], [], []
for s in args.samples:
    vcf = f"{args.results}/vcf/{s}.vcf.gz"
    callable_bed = f"{args.results}/coverage/{s}.callable.bed"
    callable_iv = {}
    for line in open(callable_bed):
        c, a, b = line.split()[:3]
        callable_iv.setdefault(c, []).append((int(a), int(b)))

    for c in callable_iv:
        callable_iv[c].sort()
    starts = {c: [x[0] for x in v] for c, v in callable_iv.items()}

    def is_callable(c, pos):  # pos is 1-based
        i = bisect.bisect_left(starts.get(c, []), pos) - 1
        return i >= 0 and callable_iv[c][i][1] >= pos

    # filtered calls, keyed by position, with the annotations the rules use
    filtered = {}
    q = "%CHROM\t%POS\t%TYPE\t%FILTER\t%QD\t%FS\t%MQ\t%MQRankSum\t%ReadPosRankSum\t%SOR\t[%GT]\n"
    for line in sh(f"bcftools query -i 'FILTER!=\"PASS\"' -f '{q}' {vcf}").splitlines():
        c, pos, typ, filt, *ann, gt = line.split("\t")
        filtered[(c, int(pos))] = (typ, filt, dict(zip(["QD", "FS", "MQ", "MQRankSum", "ReadPosRankSum", "SOR"],
                                                       map(num, ann))), gt)

    for vtype, sel in [("snv", "snps"), ("indel", "indels")]:
        d = f"{work}/{s}.exome.{vtype}"
        if not os.path.exists(f"{d}/fn.vcf.gz"):
            sh(f"rm -rf {d} && bcftools view --no-version -v {sel} {TRUTH} -Oz -o {d}.truth.vcf.gz "
               f"&& tabix -f {d}.truth.vcf.gz && bcftools view --no-version -v {sel} {vcf} -Oz -o {d}.calls.vcf.gz "
               f"&& tabix -f {d}.calls.vcf.gz && rtg vcfeval -b {d}.truth.vcf.gz -c {d}.calls.vcf.gz -t {SDF} "
               f"--bed-regions {exome} --evaluation-regions {TRUTH_BED} --squash-ploidy --threads 4 -o {d}")
        fns = [l.split("\t")[:2] for l in sh(f"bcftools query -f '%CHROM\t%POS\n' {d}/fn.vcf.gz").splitlines()]
        cls = Counter()
        rules = Counter()
        for c, pos in fns:
            pos = int(pos)
            depth = "callable" if is_callable(c, pos) else "below_min_depth"
            f = filtered.get((c, pos))
            if f:
                cls[("filtered", depth)] += 1
                if vtype == "snv":
                    fired = [r for r, test in SNP_RULES.items() if test(f[2])]
                    gt = "hom_alt" if f[3] in ("1/1", "1|1") else "het" if f[3] in ("0/1", "0|1") else "other"
                    rules[("+".join(fired) or f[1], gt, depth)] += 1
            else:
                cls[("not_called", depth)] += 1
        for (k, dep), n in sorted(cls.items()):
            fn_rows.append([s, vtype, k, dep, n, len(fns)])
        for (r, gt, dep), n in sorted(rules.items(), key=lambda x: -x[1]):
            rule_rows.append([s, r, gt, dep, n])

        if vtype == "indel":
            # longest homopolymer starting within 1 bp after the indel's anchor base
            fps = [l.split("\t") for l in sh(f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' {d}/fp.vcf.gz")
                   .splitlines()]
            regions = f"{d}/fp_regions.txt"
            with open(regions, "w") as fh:
                for c, pos, *_ in fps:
                    fh.write(f"{c}:{int(pos) + 1}-{int(pos) + 30}\n")
            seqs, cur = [], None
            for line in sh(f"samtools faidx {REF} -r {regions}").splitlines():
                if line.startswith(">"):
                    if cur is not None:
                        seqs.append(cur)
                    cur = ""
                else:
                    cur += line.upper()
            seqs.append(cur)
            bins = Counter()
            for (c, pos, ref, alt), sq in zip(fps, seqs):
                run = len(sq) - len(sq.lstrip(sq[0]))
                b = "1-3" if run < 4 else "4-7" if run < 8 else ">=8"
                kind = "del_1bp" if len(ref) - len(alt) == 1 else "ins_1bp" if len(alt) - len(ref) == 1 else "other"
                bins[(kind, b)] += 1
            for (kind, b), n in sorted(bins.items()):
                fp_rows.append([s, kind, b, n, len(fps)])

write_tsv(f"{out}/diagnosis_false_negatives{sfx}.tsv",
          ["sample", "type", "status", "depth", "n", "total_fn"], fn_rows)
write_tsv(f"{out}/diagnosis_snv_filter_rules{sfx}.tsv",
          ["sample", "rules_failed", "genotype", "depth", "n"], rule_rows)
write_tsv(f"{out}/diagnosis_indel_false_positives{sfx}.tsv",
          ["sample", "kind", "homopolymer_run", "n", "total_fp"], fp_rows)
for r in fn_rows:
    print("\t".join(map(str, r)))
