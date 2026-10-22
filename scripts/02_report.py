"""Summarise a pipeline run against config/acceptance.yaml.

Reads results/ (vcfeval summaries, panel depth, reproducibility tables) and
writes tables plus a pass/fail summary to docs/reports/<name>/.

    python scripts/02_report.py --name 01-initial-validation
"""
import argparse
import bisect
import csv
import glob
import gzip
import math
import os
import re
import shutil
from collections import defaultdict

import yaml

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
p = argparse.ArgumentParser()
p.add_argument("--name", required=True)
p.add_argument("--results", default="results")
args = p.parse_args()
os.chdir(ROOT)
acc = yaml.safe_load(open("config/acceptance.yaml"))
out = f"docs/reports/{args.name}"
os.makedirs(out, exist_ok=True)


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    ph = k / n
    d = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / d
    h = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)


def vcfeval_none_row(path):
    """The 'None' threshold row: all PASS calls, no score cut-off."""
    for line in open(path):
        f = line.split()
        if f and f[0] == "None":
            return dict(tp_base=int(f[1]), tp_call=int(f[2]), fp=int(f[3]), fn=int(f[4]))
    return dict(tp_base=0, tp_call=0, fp=0, fn=0)  # no variants in region


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


checks = []  # (criterion, sample, observed, threshold, pass)

# accuracy
acc_rows = []
for f in sorted(glob.glob(f"{args.results}/vcfeval/*.genotype.summary.txt")):
    sample, region, vtype = os.path.basename(f).split(".")[:3]
    g = vcfeval_none_row(f)
    a = vcfeval_none_row(f.replace(".genotype.", ".allele."))
    n_sens, n_prec = g["tp_base"] + g["fn"], g["tp_call"] + g["fp"]
    sens = g["tp_base"] / n_sens if n_sens else float("nan")
    prec = g["tp_call"] / n_prec if n_prec else float("nan")
    gtc = g["tp_base"] / a["tp_base"] if a["tp_base"] else float("nan")
    s_lo, s_hi = wilson(g["tp_base"], n_sens)
    p_lo, p_hi = wilson(g["tp_call"], n_prec)
    acc_rows.append([sample, region, vtype, g["tp_base"], g["fn"], g["tp_call"], g["fp"],
                     f"{sens:.4f}", f"{s_lo:.4f}", f"{s_hi:.4f}", f"{prec:.4f}", f"{p_lo:.4f}", f"{p_hi:.4f}",
                     a["tp_base"], f"{gtc:.4f}"])
    th = acc["accuracy"][vtype]
    for name, val, t in [("sensitivity", sens, th["sensitivity"]), ("precision", prec, th["precision"]),
                         ("genotype_concordance", gtc, acc["accuracy"]["genotype_concordance"])]:
        if math.isnan(val):
            why = "no truth variants" if name != "genotype_concordance" else "no allele-matched calls"
            checks.append((f"{region} {vtype} {name}", sample, why, t, "not assessable"))
        else:
            checks.append((f"{region} {vtype} {name}", sample, f"{val:.4f}", t, "PASS" if val >= t else "FAIL"))
write_tsv(f"{out}/accuracy.tsv",
          ["sample", "region", "type", "tp_baseline", "fn", "tp_call", "fp", "sensitivity", "sens_ci_low",
           "sens_ci_high", "precision", "prec_ci_low", "prec_ci_high", "tp_baseline_allele", "genotype_concordance"],
          acc_rows)

# coverage over panel targets, per gene
genes = defaultdict(list)
for line in open("data/panels/panel_targets_by_gene.bed"):
    c, s, e, g = line.split()
    genes[c].append((int(s), int(e), g))
starts = {c: [x[0] for x in sorted(v)] for c, v in genes.items()}
for c in genes:
    genes[c].sort()

min_d = acc["coverage"]["min_depth"]
cov_rows, low_rows = [], []
for f in sorted(glob.glob(f"{args.results}/coverage/*.panel_depth.tsv.gz")):
    sample = os.path.basename(f).split(".")[0]
    tot, ok = defaultdict(int), defaultdict(int)
    depths = defaultdict(list)
    run = None  # current low-coverage run: [chrom, start, end, gene, min_depth]
    with gzip.open(f, "rt") as fh:
        for line in fh:
            c, pos, d = line.split("\t")
            pos0, d = int(pos) - 1, int(d)
            i = bisect.bisect_right(starts[c], pos0) - 1
            # a base can sit in more than one gene's padded CDS; count it for each
            hits = [x[2] for x in genes[c][max(0, i - 3):i + 1] if x[0] <= pos0 < x[1]]
            for g in hits:
                tot[g] += 1
                ok[g] += d >= min_d
                depths[g].append(d)
            if d < min_d:
                g = ",".join(hits)
                if run and run[0] == c and run[2] == pos0 and run[3] == g:
                    run[2], run[4] = pos0 + 1, min(run[4], d)
                else:
                    if run:
                        low_rows.append([sample] + run)
                    run = [c, pos0, pos0 + 1, g, d]
    if run:
        low_rows.append([sample] + run)
    T, K = sum(tot.values()), sum(ok.values())
    for g in sorted(tot):
        ds = sorted(depths[g])
        cov_rows.append([sample, g, tot[g], ok[g], f"{ok[g] / tot[g]:.4f}", ds[len(ds) // 2]])
    frac = K / T
    t = acc["coverage"]["min_fraction_of_targets_at_min_depth"]
    cov_rows.append([sample, "ALL", T, K, f"{frac:.4f}", ""])
    checks.append((f"panel bases >= {min_d}x", sample, f"{frac:.4f}", t, "PASS" if frac >= t else "FAIL"))
write_tsv(f"{out}/coverage_by_gene.tsv", ["sample", "gene", "bases", f"bases_ge_{min_d}x", "fraction", "median_depth"],
          cov_rows)
write_tsv(f"{out}/low_coverage_regions.tsv", ["sample", "chrom", "start", "end", "gene", "min_depth"], low_rows)

# reproducibility
for f in glob.glob(f"{args.results}/reproducibility/within_run.*.tsv"):
    r = next(csv.DictReader(open(f), delimiter="\t"))
    same = r["differing_lines"] == "0" and r["md5_a"] == r["md5_b"]
    checks.append(("within-run identical VCF records", f"{r['run_a']} vs {r['run_b']}",
                   f"{r['differing_lines']} differing lines", "0", "PASS" if same else "FAIL"))
    shutil.copy(f, out)
rep_rows = []
for f in sorted(glob.glob(f"{args.results}/reproducibility/between_library.*.summary.txt")):
    m = re.match(r"between_library\.(.+)_vs_(.+)\.(exome|panel)\.summary\.txt", os.path.basename(f))
    a, b, region = m.groups()
    r = vcfeval_none_row(f)
    n = r["tp_base"] + r["fp"] + r["fn"]
    conc = r["tp_base"] / n if n else float("nan")
    rep_rows.append([a, b, region, r["tp_base"], r["fp"], r["fn"], f"{conc:.4f}"])
    t = acc["reproducibility"]["between_library"]["genotype_concordance"]
    if n == 0:
        checks.append((f"between-library concordance ({region})", f"{a} vs {b}", "no comparable calls", t,
                       "not assessable"))
    else:
        checks.append((f"between-library concordance ({region})", f"{a} vs {b}", f"{conc:.4f}", t,
                       "PASS" if conc >= t else "FAIL"))
write_tsv(f"{out}/between_library.tsv", ["library_a", "library_b", "region", "matched", "b_only", "a_only",
                                         "concordance"], rep_rows)

write_tsv(f"{out}/acceptance_checks.tsv", ["criterion", "sample", "observed", "threshold", "result"], checks)
for c in checks:
    print("\t".join(map(str, c)))
