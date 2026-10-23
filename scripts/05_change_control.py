"""Apply the change-control rules in config/acceptance.yaml to a pipeline change.

Needs, for the change report folder:
  test_pipeline.tsv      from 04_split_eval.py --set test --dirs <before> <after>
  acceptance_checks.tsv  from 02_report.py on the after results
and the before report's acceptance_checks.tsv.

Writes change_control.tsv (one row per rule) and panel_discordant_calls.tsv
(every PASS call that differs between the two versions in the panel targets).

    python scripts/05_change_control.py --before-report docs/reports/01-initial-validation \
        --after-report docs/reports/02-snv-filter-change --before results --after results_v2
"""
import argparse
import csv
import os
import subprocess

import yaml

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
p = argparse.ArgumentParser()
p.add_argument("--before-report", required=True)
p.add_argument("--after-report", required=True)
p.add_argument("--before", required=True)
p.add_argument("--after", required=True)
p.add_argument("--samples", nargs="+", default=["NIST7035_run1", "NIST7035_run2", "NIST7086_run1"])
args = p.parse_args()
os.chdir(ROOT)
cc = yaml.safe_load(open("config/acceptance.yaml"))["change_control"]
TRUTH = "data/truth/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED = "data/truth/HG001_GRCh37_1_22_v4.2.1_benchmark.bed"
PANEL = "data/panels/panel_targets.bed"


def sh(cmd):
    return subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True, capture_output=True, text=True).stdout


def read_tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


rules = []  # rule, detail, before, after, result

# 1. no criterion that passed before fails after
before = {(r["criterion"], r["sample"]): r for r in read_tsv(f"{args.before_report}/acceptance_checks.tsv")}
after = {(r["criterion"], r["sample"]): r for r in read_tsv(f"{args.after_report}/acceptance_checks.tsv")}
for k, b in before.items():
    a = after.get(k)
    if b["result"] == "PASS":
        ok = a is not None and a["result"] == "PASS"
        rules.append(["1 no passing criterion fails", f"{k[0]} ({k[1]})", b["observed"],
                      a["observed"] if a else "missing", "PASS" if ok else "FAIL"])

# 2. no sensitivity or precision drop above max_drop on the test set
test = read_tsv(f"{args.after_report}/test_pipeline.tsv")
by = {(r["label"], r["sample"], r["type"]): r for r in test}
for (label, s, t), r in sorted(by.items()):
    if label != args.before:
        continue
    a = by[(args.after, s, t)]
    for m in ["sensitivity", "precision"]:
        d = float(a[m]) - float(r[m])
        rules.append([f"2 {m} drop <= {cc['max_drop']} (test set)", f"{s} {t}", r[m], a[m],
                      "PASS" if d >= -cc["max_drop"] else "FAIL"])

# 4. within-run reproducibility after the change
for k, a in after.items():
    if k[0].startswith("within-run"):
        rules.append(["4 within-run identical", k[1], before[k]["observed"], a["observed"], a["result"]])

# 3. list every PASS call in the panel targets that differs between versions
disc = []
truth_gt = {}
for line in sh(f"bcftools query -R {PANEL} -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\n' {TRUTH}").splitlines():
    c, pos, ref, alt, gt = line.split("\t")
    truth_gt[(c, pos, ref, alt)] = gt
highconf = sh(f"bedtools intersect -a {PANEL} -b {TRUTH_BED} | cut -f1-3")
hc = [(l.split()[0], int(l.split()[1]), int(l.split()[2])) for l in highconf.splitlines()]
q = "%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%INFO/DP\\t%INFO/SOR\\t[%GT]\\n"
for s in args.samples:
    calls = {}
    for tag, d in [("before", args.before), ("after", args.after)]:
        for line in sh(f"bcftools query -R {PANEL} -f '{q}' {d}/vcf/{s}.vcf.gz").splitlines():
            c, pos, ref, alt, filt, dp, sor, gt = line.split("\t")
            calls.setdefault((c, pos, ref, alt), {})[tag] = (filt, dp, sor, gt)
    for k, v in sorted(calls.items(), key=lambda x: (x[0][0], int(x[0][1]))):
        fb = v.get("before", ("absent",))[0]
        fa = v.get("after", ("absent",))[0]
        if (fb == "PASS") != (fa == "PASS") or v.get("before", (0, 0, 0, None))[3] != v.get("after", (0, 0, 0, None))[3]:
            info = v.get("after") or v.get("before")
            in_hc = any(c == k[0] and a < int(k[1]) <= b for c, a, b in hc)
            truth = truth_gt.get(k, "absent from truth" if in_hc else "outside high-confidence")
            disc.append([s, *k, fb, fa, info[3], info[1], info[2], truth])
write_tsv(f"{args.after_report}/panel_discordant_calls.tsv",
          ["sample", "chrom", "pos", "ref", "alt", "filter_before", "filter_after", "gt", "dp", "sor", "truth_gt"], disc)
rules.append(["3 panel discordant calls listed", f"{len(disc)} rows in panel_discordant_calls.tsv", "", "",
              "review"])

write_tsv(f"{args.after_report}/change_control.tsv", ["rule", "detail", "before", "after", "result"], rules)
for r in rules:
    print("\t".join(map(str, r)))
fails = [r for r in rules if r[4] == "FAIL"]
print(f"\n{len(fails)} rule(s) failed")
