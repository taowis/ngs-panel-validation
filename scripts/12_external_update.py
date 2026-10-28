"""Verify an update of an externally supplied variant caller (report 05).

Works from the exported VCFs only, so the caller is treated as a black box.
Each version is given as one or more VCFs of the same sample and input: the
first is scored, the others are repeat runs used to measure run-to-run noise.

    python scripts/12_external_update.py --genome HG003 \
        --before 1.9.0=work_dv/gadi/results/deepvariant_1.9.0_run1/HG003.vcf.gz,<repeat> \
        --after 1.10.0=<run1>,<repeat> \
        --context gatk_v2=results_v2_HG003/vcf/HG003.vcf.gz \
        --out docs/reports/05-deepvariant-update

Writes accuracy.tsv (each call set, panel and exome evaluable regions, SNVs and
indels, genotype and allele mode), repeatability.tsv, change_control.tsv (the
rules in config/acceptance.yaml), panel_discordant_calls.tsv (every PASS call
in the panel targets that differs between the two versions, with the truth) and
record_differences.tsv (every record of the two VCFs, PASS or not, counted by
the fields that changed, so relabelled filters and changed scores are seen too).
"""
import argparse
import os

import yaml

from genomes import (PANEL, PANEL_BY_GENE, ROOT, exome_bed, load, rates, sh, vcfeval, wilson,
                     write_tsv)

p = argparse.ArgumentParser()
p.add_argument("--genome", required=True)
p.add_argument("--before", required=True, help="label=vcf[,repeat vcf ...]")
p.add_argument("--after", required=True)
p.add_argument("--context", nargs="*", default=[], help="label=vcf, scored but not part of the rules")
p.add_argument("--out", required=True)
p.add_argument("--work", default="work_dv/eval")
args = p.parse_args()
os.chdir(ROOT)
os.makedirs(args.out, exist_ok=True)
os.makedirs(args.work, exist_ok=True)
acc = yaml.safe_load(open("config/acceptance.yaml"))
g = load(args.genome)


def parse(s):
    label, files = s.split("=", 1)
    return label, files.split(",")


before, after = parse(args.before), parse(args.after)
context = [parse(c) for c in args.context]


def pass_only(vcf, name):
    out = f"{args.work}/{name}.pass.vcf.gz"
    if not os.path.exists(out + ".tbi"):
        sh(f"bcftools view --no-version -f PASS,. '{vcf}' -Oz -o {out} && tabix -f {out}")
    return out


regions = {"panel": g["evaluable"], "exome": exome_bed(g, f"{args.work}/exome_evaluable.{args.genome}.bed")}
scored = {}
rows = []
for label, files in [before, after] + context:
    vcf = pass_only(files[0], f"{label}_run1")
    scored[label] = vcf
    for region, bed in regions.items():
        for vtype in ("snv", "indel"):
            gt = vcfeval(g["truth_vcf"], vcf, bed, g["truth_bed"], f"{args.work}/{label}.{region}.{vtype}.gt", vtype)
            al = vcfeval(g["truth_vcf"], vcf, bed, g["truth_bed"], f"{args.work}/{label}.{region}.{vtype}.al", vtype,
                         squash=True)
            tpb, tpc, fp, fn = gt
            sens, prec = rates(*gt)
            slo, shi = wilson(tpb, tpb + fn)
            plo, phi = wilson(tpc, tpc + fp)
            conc = gt[0] / al[0] if al[0] else float("nan")
            rows.append([label, region, vtype, tpb, fp, fn, f"{sens:.4f}", f"{slo:.4f}-{shi:.4f}",
                         f"{prec:.4f}", f"{plo:.4f}-{phi:.4f}", f"{conc:.4f}"])
write_tsv(f"{args.out}/accuracy.tsv", ["call_set", "region", "type", "tp", "fp", "fn", "sensitivity", "sens_95ci",
                                        "precision", "prec_95ci", "genotype_concordance"], rows)
res = {(r[0], r[1], r[2]): r for r in rows}


def records(vcf, bed):
    """PASS records in bed as CHROM POS REF ALT GT, after splitting multiallelic sites"""
    out = sh(f"bcftools view --no-version -R {bed} '{vcf}' | bcftools norm -m- --no-version "
             f"| bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\n'")
    return set(out.splitlines())


# run-to-run noise: repeat runs of each version against its first run, over the exome evaluable region
rep_rows = []
for label, files in (before, after):
    first = records(scored[label], regions["exome"])
    for i, f in enumerate(files[1:], start=2):
        other = records(pass_only(f, f"{label}_run{i}"), regions["exome"])
        rep_rows.append([label, f"run1_vs_run{i}", len(first), len(other), len(first ^ other)])
write_tsv(f"{args.out}/repeatability.tsv", ["call_set", "comparison", "records_a", "records_b", "differing_records"],
          rep_rows)

# change-control rules
cc = acc["change_control"]
th = acc["accuracy"]
rules = []
for region in ("panel", "exome"):
    for vtype in ("snv", "indel"):
        b, a = res[(before[0], region, vtype)], res[(after[0], region, vtype)]
        for i, metric in ((6, "sensitivity"), (8, "precision")):
            limit = th[vtype][metric]
            bp, ap = float(b[i]) >= limit, float(a[i]) >= limit
            if bp:
                rules.append(["1 no passing criterion fails", f"{region} {vtype} {metric} >= {limit}", b[i], a[i],
                              "PASS" if ap else "FAIL"])
            drop = float(b[i]) - float(a[i])
            rules.append(["2 no drop above max_drop", f"{region} {vtype} {metric}", b[i], a[i],
                          "PASS" if drop <= cc["max_drop"] else "FAIL"])
for r in rep_rows:
    rules.append(["3 repeat runs identical", f"{r[0]} {r[1]}", r[2], r[3], "PASS" if r[4] == 0 else f"FAIL ({r[4]})"])

# every PASS call in the panel targets that differs between the versions, with the truth genotype
genes = sh(f"cat {PANEL_BY_GENE}").splitlines()


def gene_at(chrom, pos):
    for line in genes:
        c, s, e, name = line.split("\t")[:4]
        if c == chrom and int(s) < pos <= int(e):
            return name
    return ""


b_rec = {tuple(r.split("\t")[:4]): r.split("\t")[4] for r in records(scored[before[0]], PANEL)}
a_rec = {tuple(r.split("\t")[:4]): r.split("\t")[4] for r in records(scored[after[0]], PANEL)}
disc = []
for k in sorted(set(b_rec) | set(a_rec), key=lambda k: (k[0].zfill(2), int(k[1]))):
    if b_rec.get(k) == a_rec.get(k):
        continue
    truth = sh(f"bcftools view --no-version -H -r {k[0]}:{k[1]} {g['truth_vcf']} | bcftools norm -m- --no-version "
               f"| bcftools query -f '%REF\\t%ALT\\t[%GT]\\n' 2>/dev/null || true")
    tgt = next((t.split("\t")[2] for t in truth.splitlines() if tuple(t.split("\t")[:2]) == k[2:]), "")
    in_conf = sh(f"echo -e '{k[0]}\\t{int(k[1]) - 1}\\t{k[1]}' | bedtools intersect -u -a - -b {g['truth_bed']}")
    disc.append([*k, gene_at(k[0], int(k[1])), b_rec.get(k, "not called"), a_rec.get(k, "not called"),
                 tgt or ("absent" if in_conf else "outside high-confidence")])
write_tsv(f"{args.out}/panel_discordant_calls.tsv",
          ["chrom", "pos", "ref", "alt", "gene", f"gt_{before[0]}", f"gt_{after[0]}", "gt_truth"], disc)
rules.append(["4 panel calls that differ listed", "PASS calls in panel targets", len(b_rec), len(a_rec),
              f"{len(disc)} listed"])
# every record of the two first runs, PASS or not, to catch changes the accuracy figures cannot see
# (relabelled filters, changed QUAL or likelihoods, records added or dropped)
FIELDS = ["QUAL", "FILTER", "GT", "other FORMAT"]


def all_records(vcf):
    rec = {}
    for line in sh(f"bcftools view --no-version -H '{vcf}'").splitlines():
        f = line.split("\t")
        fmt = dict(zip(f[8].split(":"), f[9].split(":")))
        gt = fmt.pop("GT", ".")
        rec[tuple(f[:2] + f[3:5])] = (f[5], f[6], gt, ":".join(f"{k}={v}" for k, v in sorted(fmt.items())))
    return rec


ra, rb = all_records(before[1][0]), all_records(after[1][0])
counts = {}
for k in set(ra) | set(rb):
    a, b = ra.get(k), rb.get(k)
    if a == b:
        continue
    if a is None or b is None:
        key = ("record only in " + (after[0] if a is None else before[0]), "", "")
    else:
        changed = ",".join(n for n, x, y in zip(FIELDS, a, b) if x != y)
        key = (changed, a[1], b[1])
    counts[key] = counts.get(key, 0) + 1
write_tsv(f"{args.out}/record_differences.tsv",
          ["fields_changed", f"filter_{before[0]}", f"filter_{after[0]}", "records"],
          [list(k) + [n] for k, n in sorted(counts.items(), key=lambda x: -x[1])])
rules.append(["5 all record differences listed", f"records in {before[0]} / {after[0]}", len(ra), len(rb),
              f"{sum(counts.values())} differ, see record_differences.tsv"])
write_tsv(f"{args.out}/change_control.tsv", ["rule", "detail", before[0], after[0], "result"], rules)
for r in rules:
    print("\t".join(map(str, r)))
