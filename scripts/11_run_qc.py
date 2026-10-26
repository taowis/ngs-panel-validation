"""Run-level QC metrics (NPAAC S6.6, S6.7) for each sample, and a check
against the acceptable ranges in config/run_qc.yaml.

reads     samtools flagstat on the duplicate-marked BAM: primary reads, % mapped,
          % properly paired; duplicate rate from MarkDuplicates; mean depth over
          the capture from mosdepth
variants  PASS calls inside the merged capture BED, multiallelic records split:
          SNV and indel counts, SNV Ti/Tv, SNV het/hom-alt ratio, and the
          fraction of SNVs and indels in dbSNP 138 (same position and alleles)

    python scripts/11_run_qc.py                      # all validation samples
    python scripts/11_run_qc.py --genomes HG002 --out /tmp/qc.tsv
    python scripts/11_run_qc.py --recheck            # ranges only, on the existing table
"""
import argparse
import os

import yaml

from genomes import GENOMES, REF, ROOT, load, read_tsv, sh, write_tsv

p = argparse.ArgumentParser()
p.add_argument("--genomes", nargs="+", default=GENOMES)
p.add_argument("--ranges", default="config/run_qc.yaml")
p.add_argument("--dbsnp", default="data/reference/dbsnp_138.b37.vcf.gz")
p.add_argument("--out", default="docs/reports/04-run-qc/run_qc.tsv")
p.add_argument("--recheck", action="store_true", help="apply the ranges to --out without recomputing")
args = p.parse_args()
os.chdir(ROOT)
ranges = yaml.safe_load(open(args.ranges)) if os.path.exists(args.ranges) else None

COLS = ["primary_reads", "pct_mapped", "pct_properly_paired", "pct_duplicates", "mean_capture_depth",
        "pass_snvs", "pass_indels", "snv_titv", "snv_het_hom", "pct_snvs_in_dbsnp", "pct_indels_in_dbsnp"]


def read_metrics(g, s):
    fs = {}
    for line in sh(f"samtools flagstat -@ 4 -O tsv {g['outdir']}/bam/{s}.md.bam").splitlines():
        n, _, name = line.split("\t")
        if not n.endswith("%"):
            fs[name] = int(n)
    with open(f"{g['outdir']}/bam/{s}.markdup_metrics.txt") as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if l.strip() and not l.startswith("#")]
    dup = float(dict(zip(lines[0], lines[1]))["PERCENT_DUPLICATION"])
    depth = next(l.split()[3] for l in open(f"{g['outdir']}/coverage/{s}.mosdepth.summary.txt")
                 if l.startswith("total_region"))
    return [fs["primary"], 100 * fs["primary mapped"] / fs["primary"],
            100 * fs["properly paired"] / fs["paired in sequencing"], 100 * dup, float(depth)]


def variant_metrics(g, s, capture):
    vcf = f"{g['outdir']}/vcf/{s}.vcf.gz"
    os.makedirs(f"{g['outdir']}/run_qc", exist_ok=True)
    tmp = f"{g['outdir']}/run_qc/{s}.pass.vcf.gz"  # annotate -a needs an indexed input
    sh(f"bcftools view -f PASS -R {capture} {vcf} | bcftools norm --no-version -m- -f {REF} -Oz -o {tmp} 2>/dev/null "
       f"&& tabix -f {tmp}")
    rows = sh(f"bcftools annotate -x ID -a {args.dbsnp} -c ID {tmp} "
              f"| bcftools query -f '%TYPE\\t%REF\\t%ALT\\t%ID[\\t%GT]\\n'")
    ts = tv = het = hom = snv_db = ind = ind_db = 0
    tsset = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    for line in rows.splitlines():
        typ, ref, alt, vid, gt = line.split("\t")
        known = vid != "."
        if typ == "SNP":
            ts, tv = (ts + 1, tv) if (ref, alt) in tsset else (ts, tv + 1)
            snv_db += known
            a = gt.replace("|", "/").split("/")
            het, hom = (het + 1, hom) if len(set(a)) > 1 else (het, hom + 1)
        elif typ == "INDEL":
            ind += 1
            ind_db += known
    snvs = ts + tv
    return [snvs, ind, ts / tv, het / hom, 100 * snv_db / snvs, 100 * ind_db / ind]


def check(assay, vals):
    if not ranges:
        return assay, "not checked"
    a = ranges["assays"][assay]
    out = []
    for k, v in zip(COLS, vals):
        lo, hi = a["ranges"].get(k, [None, None])
        if (lo is not None and v < lo) or (hi is not None and v > hi):
            out.append(k)
    return a["name"], ",".join(out) or "none"


rows = []
if args.recheck:
    names = {a["name"]: k for k, a in ranges["assays"].items()}
    for r in read_tsv(args.out):
        name, out = check(names.get(r["assay"], r["assay"]), [float(r[k]) for k in COLS])
        rows.append([r["genome"], r["sample"], name] + [r[k] for k in COLS] + [out])
for gname in [] if args.recheck else args.genomes:
    g = load(gname)
    assay = os.path.basename(g["capture"])
    capture = f"{g['outdir']}/coverage/capture.merged.bed"
    if not os.path.exists(capture):
        sh(f"gzip -dc {g['capture']} | sed 's/^chr//' | grep -v '^track\\|^browser' | cut -f1-3 "
           f"| sort -k1,1 -k2,2n | bedtools merge > {capture}")
    for s in g["samples"]:
        if not os.path.exists(f"{g['outdir']}/bam/{s}.md.bam"):
            continue  # NIST7035_run2: BAM not kept, same FASTQ as run1
        vals = read_metrics(g, s) + variant_metrics(g, s, capture)
        name, out = check(assay, vals)
        rows.append([gname, s, name] + [f"{v:.4g}" if isinstance(v, float) else v for v in vals] + [out])
for row in rows:
    print("\t".join(map(str, row)), flush=True)
write_tsv(args.out, ["genome", "sample", "assay"] + COLS + ["out_of_range"], rows)
