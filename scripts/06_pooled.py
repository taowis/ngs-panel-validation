"""Report 03, analyses 2 and 5: panel accuracy pooled over the five genomes, and
accuracy against depth for the downsampled HG002.

Reads accuracy.tsv and coverage_by_gene.tsv written by 02_report.py into the
report 03 folder (HG002 to HG005, and the HG002 depth levels) and report 02's
accuracy.tsv for HG001.

    python scripts/06_pooled.py
"""
import argparse
import os
import re

from genomes import GENOMES, ROOT, load, read_tsv, sh, wilson, write_tsv

p = argparse.ArgumentParser()
p.add_argument("--report", default="docs/reports/03-more-reference-samples")
p.add_argument("--hg001-report", default="docs/reports/02-snv-filter-change")
args = p.parse_args()
os.chdir(ROOT)

acc = read_tsv(f"{args.report}/accuracy.tsv") + read_tsv(f"{args.hg001_report}/accuracy.tsv")
by = {(r["sample"], r["region"], r["type"]): r for r in acc}

# 2. pooled panel accuracy, one sample per genome
rows, per_genome = [], []
for vtype, sel in [("snv", "snps"), ("indel", "indels")]:
    tot = dict(tp_baseline=0, fn=0, tp_call=0, fp=0)
    distinct = set()
    for gname in GENOMES:
        g = load(gname)
        r = by[(g["main_sample"], "panel", vtype)]
        for k in tot:
            tot[k] += int(r[k])
        per_genome.append([gname, g["main_sample"], vtype, r["tp_baseline"], r["fn"], r["tp_call"], r["fp"]])
        q = sh(f"bcftools view --no-version -H -v {sel} -R {g['evaluable']} {g['truth_vcf']} | cut -f1,2,4,5 | tr '\\t' :")
        distinct |= set(q.split())
    n_sens, n_prec = tot["tp_baseline"] + tot["fn"], tot["tp_call"] + tot["fp"]
    sens = tot["tp_baseline"] / n_sens if n_sens else float("nan")
    prec = tot["tp_call"] / n_prec if n_prec else float("nan")
    s_lo, s_hi = wilson(tot["tp_baseline"], n_sens)
    p_lo, p_hi = wilson(tot["tp_call"], n_prec)
    rows.append([vtype, len(GENOMES), tot["tp_baseline"], tot["fn"], tot["tp_call"], tot["fp"], f"{sens:.4f}",
                 f"{s_lo:.4f}", f"{s_hi:.4f}", f"{prec:.4f}", f"{p_lo:.4f}", f"{p_hi:.4f}", n_sens, len(distinct)])
write_tsv(f"{args.report}/pooled_panel_accuracy.tsv",
          ["type", "genomes", "tp_baseline", "fn", "tp_call", "fp", "sensitivity", "sens_ci_low", "sens_ci_high",
           "precision", "prec_ci_low", "prec_ci_high", "truth_variants", "distinct_truth_variants"], rows)
write_tsv(f"{args.report}/pooled_panel_inputs.tsv",
          ["genome", "sample", "type", "tp_baseline", "fn", "tp_call", "fp"], per_genome)

# 5. accuracy against depth, HG002
g = load("HG002")
cov = {r["sample"]: r for r in read_tsv(f"{args.report}/coverage_by_gene.tsv") if r["gene"] == "ALL"}
levels = sorted({r["sample"] for r in acc if re.fullmatch(r"HG002(_\d+x)?", r["sample"])},
                key=lambda s: int(s.split("_")[1][:-1]) if "_" in s else 10 ** 6)
drows = []
for s in levels:
    mean = None
    for line in open(f"{g['outdir']}/coverage/{s}.mosdepth.summary.txt"):
        f = line.split()
        if f[0] == "total_region":
            mean = float(f[3])
    for region in ["exome", "panel"]:
        for vtype in ["snv", "indel"]:
            r = by[(s, region, vtype)]
            drows.append([s, f"{mean:.1f}", cov[s]["fraction"], region, vtype, r["tp_baseline"], r["fn"],
                          r["sensitivity"], r["sens_ci_low"], r["precision"]])
write_tsv(f"{args.report}/depth_downsampling.tsv",
          ["sample", "mean_depth_capture", "panel_fraction_ge_20x", "region", "type", "tp_baseline", "fn",
           "sensitivity", "sens_ci_low", "precision"], drows)
for r in rows + drows:
    print("\t".join(map(str, r)))
