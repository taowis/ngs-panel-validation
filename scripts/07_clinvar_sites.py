"""Report 03, analysis 4: depth at ClinVar pathogenic and likely pathogenic sites
in the panel targets.

Small variants only (REF and ALT up to 50 bp, the scope of SOP-01); longer
records are counted and left out. A site counts as covered when every reference base it spans inside the targets
is at >= min_depth (samtools depth from the pipeline: MAPQ >= 20, base quality
>= 20, duplicates excluded). "In capture" means inside the kit's capture BED.

    python scripts/07_clinvar_sites.py --clinvar data/reference/clinvar_20261019.vcf.gz
"""
import argparse
import bisect
import gzip
import os
from collections import defaultdict

import yaml

from genomes import GENOMES, PANEL, PANEL_BY_GENE, ROOT, load, sh, write_tsv

p = argparse.ArgumentParser()
p.add_argument("--clinvar", required=True)
p.add_argument("--report", default="docs/reports/03-more-reference-samples")
p.add_argument("--genomes", nargs="+", default=GENOMES)
args = p.parse_args()
os.chdir(ROOT)
min_d = yaml.safe_load(open("config/acceptance.yaml"))["coverage"]["min_depth"]
PLP = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}


def intervals(lines):
    iv = defaultdict(list)
    for line in lines:
        f = line.split("\t")
        iv[f[0].removeprefix("chr")].append((int(f[1]), int(f[2]), f[3].strip() if len(f) > 3 else ""))
    for c in iv:
        iv[c].sort()
    return iv, {c: [x[0] for x in v] for c, v in iv.items()}


def hits(iv, starts, c, pos0):
    i = bisect.bisect_right(starts.get(c, []), pos0)
    return [x for x in iv.get(c, [])[max(0, i - 4):i] if x[0] <= pos0 < x[1]]


genes, gstarts = intervals(open(PANEL_BY_GENE))
targets, tstarts = intervals(open(PANEL))

# sites: positions are 0-based; span = reference bases inside the targets
sites, n_long = [], 0
q = "%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\t%INFO/CLNSIG\\n"
for line in sh(f"bcftools query -R {PANEL} -f '{q}' {args.clinvar}").splitlines():
    c, pos, ref, alt, vid, sig = line.split("\t")
    if sig.split("|")[0] not in PLP:
        continue
    if max(len(ref), *(len(a) for a in alt.split(","))) > 50:
        n_long += 1
        continue
    p0 = int(pos) - 1
    span = [x for x in range(p0, p0 + len(ref)) if hits(targets, tstarts, c, x)]
    if not span:
        continue
    gene = ",".join(sorted({x[2] for b in span for x in hits(genes, gstarts, c, b)}))
    sites.append((c, int(pos), ref, alt, vid, sig, gene, span))
need = defaultdict(set)
for s in sites:
    need[s[0]].update(s[7])

summary, by_gene, low = [], [], []
for gname in args.genomes:
    g = load(gname)
    cap, cstarts = intervals(gzip.open(g["capture"], "rt") if g["capture"].endswith(".gz") else open(g["capture"]))
    for sample in g["samples"]:
        if sample.endswith("_run2"):  # same data as run 1
            continue
        depth = {}
        with gzip.open(f"{g['outdir']}/coverage/{sample}.panel_depth.tsv.gz", "rt") as fh:
            for line in fh:
                c, pos, d = line.split("\t")
                if int(pos) - 1 in need[c]:
                    depth[(c, int(pos) - 1)] = int(d)
        n_ok = n_cap = 0
        gene_n, gene_ok = defaultdict(int), defaultdict(int)
        for c, pos, ref, alt, vid, sig, gene, span in sites:
            dmin = min(depth.get((c, b), 0) for b in span)
            in_cap = any(hits(cap, cstarts, c, b) for b in span)
            ok = dmin >= min_d
            n_ok += ok
            n_cap += in_cap
            gene_n[gene] += 1
            gene_ok[gene] += ok
            if not ok:
                low.append([gname, sample, gene, c, pos, ref, alt, vid, sig, dmin, "yes" if in_cap else "no"])
        summary.append([gname, sample, len(sites), n_cap, n_ok, f"{n_ok / len(sites):.4f}"])
        for gene in sorted(gene_n):
            by_gene.append([gname, sample, gene, gene_n[gene], gene_ok[gene], f"{gene_ok[gene] / gene_n[gene]:.4f}"])

write_tsv(f"{args.report}/clinvar_sites_summary.tsv",
          ["genome", "sample", "plp_sites", "in_capture", f"at_ge_{min_d}x", "fraction", "longer_than_50bp_excluded"],
          [r + [n_long] for r in summary])
write_tsv(f"{args.report}/clinvar_sites_by_gene.tsv",
          ["genome", "sample", "gene", "plp_sites", f"at_ge_{min_d}x", "fraction"], by_gene)
write_tsv(f"{args.report}/clinvar_sites_below_{min_d}x.tsv",
          ["genome", "sample", "gene", "chrom", "pos", "ref", "alt", "clinvar_id", "clnsig", "min_depth",
           "in_capture"], low)
print(f"{len(sites)} small P/LP sites in the targets; {n_long} records longer than 50 bp left out")
for r in summary:
    print("\t".join(map(str, r)))
