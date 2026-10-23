"""Report 03, analysis 6: exome-wide accuracy by GIAB difficult-region
stratification, and PMS2 (pseudogene PMS2CL) against MLH1 for comparison.

One sample per genome. Strata are GIAB v3.6 GRCh37 BEDs in
data/reference/stratifications/, intersected with that genome's exome
evaluable region.

    python scripts/08_strata.py
"""
import argparse
import os

from genomes import GENOMES, PANEL_BY_GENE, ROOT, exome_bed, load, rates, read_tsv, sh, vcfeval, write_tsv

p = argparse.ArgumentParser()
p.add_argument("--report", default="docs/reports/03-more-reference-samples")
p.add_argument("--genomes", nargs="+", default=GENOMES)
p.add_argument("--hg001-report", default="docs/reports/02-snv-filter-change")
args = p.parse_args()
os.chdir(ROOT)
D = "data/reference/stratifications"
STRATA = {
    "not_difficult": f"{D}/GRCh37_notinalldifficultregions.bed.gz",
    "all_difficult": f"{D}/GRCh37_alldifficultregions.bed.gz",
    "homopolymer": f"{D}/GRCh37_AllHomopolymers_ge7bp_imperfectge11bp_slop5.bed.gz",
    "segdup": f"{D}/GRCh37_segdups.bed.gz",
    "lowmap_or_segdup": f"{D}/GRCh37_alllowmapandsegdupregions.bed.gz",
    "gc_lt25_or_gt65": f"{D}/GRCh37_gclt25orgt65_slop50.bed.gz",
}

rows = []
for gname in args.genomes:
    g = load(gname)
    s = g["main_sample"]
    work = f"{g['outdir']}/strata"
    os.makedirs(work, exist_ok=True)
    exome = exome_bed(g, f"{work}/exome_evaluable.bed")
    for name, bed in STRATA.items():
        region = f"{work}/{name}.bed"
        if not os.path.exists(region):
            sh(f"gzip -dc {bed} | grep -v '^#' | sed 's/^chr//' | cut -f1-3 "
               f"| bedtools intersect -a {exome} -b - | sort -k1,1 -k2,2n | bedtools merge > {region}")
        bp = sum(int(e) - int(b) for _, b, e in (l.split()[:3] for l in open(region)))
        for vtype in ["snv", "indel"]:
            tpb, tpc, fp, fn = vcfeval(g["truth_vcf"], f"{g['outdir']}/vcf/{s}.vcf.gz", region, g["truth_bed"],
                                       f"{work}/{s}.{name}.{vtype}", vtype)
            sens, prec = rates(tpb, tpc, fp, fn)
            rows.append([gname, s, name, bp, vtype, tpb, fn, tpc, fp, f"{sens:.4f}", f"{prec:.4f}"])
write_tsv(f"{args.report}/strata_accuracy.tsv",
          ["genome", "sample", "stratum", "evaluable_bp", "type", "tp_baseline", "fn", "tp_call", "fp",
           "sensitivity", "precision"], rows)

# PMS2 and MLH1: coding sequence in high-confidence, bases at >= 20x, reads with MAPQ < 20
cov = {(r["sample"], r["gene"]): r for f in [args.report, args.hg001_report]
       for r in read_tsv(f"{f}/coverage_by_gene.tsv")}
grows = []
for gname in args.genomes:
    g = load(gname)
    s = g["main_sample"]
    work = f"{g['outdir']}/strata"
    for gene in ["PMS2", "MLH1"]:
        bed = f"{work}/{gene}.bed"
        sh(f"awk '$4 == \"{gene}\"' {PANEL_BY_GENE} | cut -f1-3 > {bed}")
        bp = sum(int(e) - int(b) for _, b, e in (l.split() for l in open(bed)))
        hc = sh(f"bedtools intersect -a {bed} -b {g['truth_bed']} | awk '{{n += $3 - $2}} END {{print n + 0}}'")
        bam = f"{g['outdir']}/bam/{s}.md.bam"
        # primary, mapped, not duplicate or QC fail
        n_all = int(sh(f"samtools view -c -F 0xF04 -L {bed} {bam}"))
        n_q20 = int(sh(f"samtools view -c -F 0xF04 -q 20 -L {bed} {bam}"))
        c = cov[(s, gene)]
        grows.append([gname, s, gene, bp, int(hc), f"{int(hc) / bp:.4f}", c["fraction"], c["median_depth"],
                      n_all, f"{1 - n_q20 / n_all:.4f}" if n_all else "nan"])
write_tsv(f"{args.report}/pms2.tsv",
          ["genome", "sample", "gene", "target_bp", "highconf_bp", "highconf_fraction", "fraction_ge_20x",
           "median_depth", "reads", "fraction_mapq_lt_20"], grows)
for r in rows + grows:
    print("\t".join(map(str, r)))
