"""Build the hereditary cancer panel BED files (GRCh37).

Genes: green (diagnostic grade) genes from the PanelApp Australia snapshots in
data/panels/panelapp/. Regions: CDS of all protein-coding Ensembl 87
transcripts, padded and merged. The evaluable region is the part of the panel
that is both captured by the exome kit and inside the GIAB high-confidence
regions; only there can sensitivity and precision be measured.

Outputs in data/panels/:
  panel_genes.tsv          gene list with source panels
  panel_targets.bed        padded CDS, merged
  panel_evaluable.bed      targets & capture & GIAB high-confidence
  panel_summary.tsv        per gene bases and truth variant counts
"""
import argparse
import glob
import gzip
import json
import os
from collections import defaultdict

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
p = argparse.ArgumentParser()
p.add_argument("--pad", type=int, default=20, help="bp added each side of every CDS block")
p.add_argument("--gtf", default="data/reference/Homo_sapiens.GRCh37.87.gtf.gz")
p.add_argument("--capture", default="data/raw/nexterarapidcapture_expandedexome_targetedregions.bed.gz")
p.add_argument("--truth-bed", default="data/truth/HG001_GRCh37_1_22_v4.2.1_benchmark.bed")
p.add_argument("--truth-vcf", default="data/truth/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz")
args = p.parse_args()
os.chdir(ROOT)
out = "data/panels"


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def merge(iv):
    """iv: dict chrom -> list of (start, end), half-open. Returns merged copy."""
    res = {}
    for c, xs in iv.items():
        xs = sorted(xs)
        m = [list(xs[0])]
        for s, e in xs[1:]:
            if s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        res[c] = [tuple(x) for x in m]
    return res


def intersect(a, b):
    res = {}
    for c in a:
        xs, ys, i, j, r = a[c], b.get(c, []), 0, 0, []
        while i < len(xs) and j < len(ys):
            s, e = max(xs[i][0], ys[j][0]), min(xs[i][1], ys[j][1])
            if s < e:
                r.append((s, e))
            if xs[i][1] < ys[j][1]:
                i += 1
            else:
                j += 1
        res[c] = r
    return res


def read_bed(path):
    iv = defaultdict(list)
    with opener(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            c, s, e = line.split("\t")[:3]
            iv[c.removeprefix("chr")].append((int(s), int(e)))
    return merge(iv)


def size(iv):
    return sum(e - s for xs in iv.values() for s, e in xs)


def write_bed(path, iv):
    with open(path, "w") as fh:
        for c in sorted(iv, key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x)):
            for s, e in iv[c]:
                fh.write(f"{c}\t{s}\t{e}\n")


# genes
genes = {}
for f in sorted(glob.glob(f"{out}/panelapp/panel_*.json")):
    d = json.load(open(f))
    for g in d["genes"]:
        if g["confidence_level"] != "3":
            continue
        gd = g["gene_data"]
        eid = gd["ensembl_genes"]["GRch37"]["87"]["ensembl_id"]
        genes.setdefault(eid, {"symbol": gd["gene_symbol"], "panels": []})
        genes[eid]["panels"].append(f"{d['name']} v{d['version']}")

# CDS per gene, half-open BED coordinates
cds = defaultdict(lambda: defaultdict(list))
with gzip.open(args.gtf, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if f[2] != "CDS" or 'transcript_biotype "protein_coding"' not in f[8]:
            continue
        eid = f[8].split('gene_id "')[1].split('"')[0]
        if eid in genes:
            cds[eid][f[0]].append((int(f[3]) - 1 - args.pad, int(f[4]) + args.pad))

missing = [genes[e]["symbol"] for e in genes if e not in cds]
if missing:
    raise SystemExit(f"no CDS found for {missing}")

capture = read_bed(args.capture)
highconf = read_bed(args.truth_bed)

# truth variants by position
truth = defaultdict(list)
with gzip.open(args.truth_vcf, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        c, pos, _, ref, alt = line.split("\t")[:5]
        kind = "snv" if len(ref) == 1 and all(len(a) == 1 for a in alt.split(",")) else "indel"
        truth[c].append((int(pos) - 1, kind))


def count(iv):
    n = {"snv": 0, "indel": 0}
    for c, xs in iv.items():
        for p0, kind in truth.get(c, []):
            if any(s <= p0 < e for s, e in xs):
                n[kind] += 1
    return n


all_t, all_e = defaultdict(list), defaultdict(list)
rows = []
for eid, g in sorted(genes.items(), key=lambda x: x[1]["symbol"]):
    t = merge(cds[eid])
    ev = intersect(intersect(t, capture), highconf)
    n = count(ev)
    rows.append([g["symbol"], eid, size(t), size(intersect(t, capture)), size(intersect(t, highconf)),
                 size(ev), n["snv"], n["indel"]])
    for c in t:
        all_t[c] += t[c]
        all_e[c] += ev.get(c, [])

all_t, all_e = merge(all_t), merge({c: x for c, x in all_e.items() if x})
write_bed(f"{out}/panel_targets.bed", all_t)
write_bed(f"{out}/panel_evaluable.bed", all_e)
n = count(all_e)
rows.append(["TOTAL", "", size(all_t), size(intersect(all_t, capture)), size(intersect(all_t, highconf)),
             size(all_e), n["snv"], n["indel"]])

with open(f"{out}/panel_genes.tsv", "w") as fh:
    fh.write("gene\tensembl_id_grch37\tsource_panels\n")
    for eid, g in sorted(genes.items(), key=lambda x: x[1]["symbol"]):
        fh.write(f"{g['symbol']}\t{eid}\t{'; '.join(g['panels'])}\n")

with open(f"{out}/panel_summary.tsv", "w") as fh:
    fh.write("gene\tensembl_id\ttarget_bp\tcaptured_bp\thighconf_bp\tevaluable_bp\ttruth_snv\ttruth_indel\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")

for r in rows:
    print("\t".join(map(str, r)))
