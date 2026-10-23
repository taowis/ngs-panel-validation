"""Report 03, analysis 8: this pipeline against the sequencing laboratory's own
calls on the same reads (Oslo University Hospital, GATK HaplotypeCaller with
VQSR, published by GIAB).

Both call sets are scored against the truth over the same exome evaluable
region (capture & high-confidence), PASS records only, and compared with each
other over that region.

    python scripts/10_interlab.py
"""
import argparse
import os

from genomes import ROOT, exome_bed, load, rates, sh, vcfeval, write_tsv

p = argparse.ArgumentParser()
p.add_argument("--report", default="docs/reports/03-more-reference-samples")
args = p.parse_args()
os.chdir(ROOT)
OSLO = {"HG002": "data/raw/oslo/HG002-HG003-HG004.jointVC.filter.vcf",
        "HG003": "data/raw/oslo/HG002-HG003-HG004.jointVC.filter.vcf",
        "HG004": "data/raw/oslo/HG002-HG003-HG004.jointVC.filter.vcf",
        "HG005": "data/raw/oslo/HG005.singleSampleVC.filter.vcf"}

acc_rows, conc_rows = [], []
for gname, src in OSLO.items():
    g = load(gname)
    work = f"{g['outdir']}/interlab"
    os.makedirs(work, exist_ok=True)
    exome = exome_bed(g, f"{work}/exome_evaluable.bed")
    names = sh(f"bcftools query -l {src}").split()
    col = [n for n in names if gname in n]
    assert len(col) == 1, f"no single {gname} column in {src}: {names}"
    theirs = f"{work}/oslo.{gname}.vcf.gz"
    if not os.path.exists(theirs + ".tbi"):
        # this sample's non-reference genotypes from the joint call set, otherwise as published
        sh(f"bcftools view --no-version -s {col[0]} -c 1 {src} | bcftools sort -Oz -o {theirs} && tabix -f {theirs}")
    ours = f"{g['outdir']}/vcf/{g['main_sample']}.vcf.gz"
    for label, vcf in [("this pipeline (v2)", ours), ("Oslo", theirs)]:
        for vtype in ["snv", "indel"]:
            tag = "ours" if label.startswith("this") else "oslo"
            tpb, tpc, fp, fn = vcfeval(g["truth_vcf"], vcf, exome, g["truth_bed"], f"{work}/{tag}.{vtype}", vtype)
            sens, prec = rates(tpb, tpc, fp, fn)
            acc_rows.append([gname, label, vtype, tpb, fn, tpc, fp, f"{sens:.4f}", f"{prec:.4f}"])
    for vtype in ["snv", "indel"]:
        # baseline = Oslo, calls = ours: tp both, fp ours only, fn Oslo only
        tpb, tpc, fp, fn = vcfeval(theirs, ours, exome, None, f"{work}/ours_vs_oslo.{vtype}", vtype)
        conc_rows.append([gname, vtype, tpb, fp, fn, f"{tpb / (tpb + fp + fn):.4f}"])
write_tsv(f"{args.report}/interlab_accuracy.tsv",
          ["genome", "call_set", "type", "tp_baseline", "fn", "tp_call", "fp", "sensitivity", "precision"], acc_rows)
write_tsv(f"{args.report}/interlab_concordance.tsv",
          ["genome", "type", "both", "this_pipeline_only", "oslo_only", "concordance"], conc_rows)
for r in acc_rows + conc_rows:
    print("\t".join(map(str, r)))
